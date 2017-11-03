#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "omp.h"

/*
 * C-function to inpaint missing data wedges in cone-beam sinogram (currently accepts 2D slices only)
 *
 * Inputs:
 * 1. 2D sinogram (A) where x-axis are detectors and y-axis are angles
 * 2. sinoMultip: a multiplier for the sinogram threshold in the range (0,1]; default: 0.3
 * 3. gradMultip: a multiplier for the gradient threshold in the range (0,1]; default: 0.05
 * 4. pad_lookup: an integer number of pixels to add to calculate a median of usabale values; default: 10
 * 5. sino_roll: an integer number to roll the sinogram from the top to bottom symmetrically; default: 200
 *
 * Outputs:
 * 1. Inpainted sinogram (B)
 * 2. dy^2 (squared gradient in Y direction) (U)
 *
 * to compile with OMP support: mex SinoConePad.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
 * gcc -shared -Wl,-soname,SinoPad -fopenmp -o SinoPad.so -fPIC SinoConePad.c
 *
 * Harwell/18.04.17
 */

float pad_func_joint(float *A, float *B, float *gradYs, int dimX,  int dimY, float edgeval, float gradval, int pad_lookup, int mid_val);
float gradYsq_func(float *U, float *gradYs, int dimX,  int dimY);
float copyAr(float *B, float *U, int dimX,  int dimY);
float copyAr_roll(float *B, float *U, int dimX,  int dimY, int roll_value, int switcher);
float cleaning_outliers_func(float *A, float *B, int dimX, int dimY, float edgeval);

/*TODO: reading parameters from Python, would be great to add some fool-proof system */
/*TODO: Passing printf messages to the user from python? */

/*Handling Python input data*/

/* A: singoram [detectors, angles] */
/* sinoMultip: a multiplier for the sinogram in the range (0,1];   */
/* gradMultip: a multiplier for the gradient in the range (0,1];   */
/* pad_lookup:  number of pixels to add to calculate a median of usabale values */
/* roll_value: symmetric sinogram roll from the top to bottom, 0 - do not roll */

void SinoPad(float *A, float sinoMultip, float gradMultip, int pad_lookup, int roll_value, int dimX, int dimY, float *B, float *gradYs) 
{
    int i, j, iter, counterElements, iterations_number, mid_val;
    float *U, edgeval, gradval, sumW;
              
    mid_val = (int)(0.5f*pad_lookup);
    
    /*TODO: some general checks are required as were in the MEX file bellow! */
    
    /* if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) {mexErrMsgTxt("The input in a single precision is required"); }
    /*if(nrhs != 5) mexErrMsgTxt("Five inputs is reqired ");    
    if ((pad_lookup <= 3) || (gradMultip > 200)) {mexErrMsgTxt("Too small or too large value for padding, try 10 "); }
    if ((roll_value < 0)) {mexErrMsgTxt(" Symmetric sinogram rolling value must be positive "); }        
    if (number_of_dims == 3) {  mexErrMsgTxt("Currently only 2D version is available"); }
    */
     if ((sinoMultip <= 0) || (sinoMultip > 1)) {
     fprintf(stderr, "The sinogram multiplier must be in the range (0,1], e.g. try: 0.3 \n");
     exit(EXIT_FAILURE); }
     if ((gradMultip <= 0) || (gradMultip > 1)) {
     fprintf(stderr, "The gradient multiplier must be in the range (0,1], e.g. try 0.05  \n");
     exit(EXIT_FAILURE); }

     U = (float*) calloc (dimY*dimX,sizeof(float));          
     /*gradYs = (float*)mxGetPr(plhs[1] = mxCreateNumericArray(2, dim_array, mxSINGLE_CLASS, mxREAL));*/
        
        /* find a critical edge-value of usable sinogram (tweakable by changing sinoMultip) */
        edgeval = 0.0f;
        for(i=0; i<dimY*dimX; i++) if (A[i] > edgeval) edgeval = A[i];
        edgeval = edgeval*sinoMultip; /* max value in sinogram multipled by sinoMultip*/

        if (roll_value > 0) copyAr_roll(A, U, dimX, dimY, roll_value, 0); /* copying A to U with symmetric rolling */
        else copyAr(A, U, dimX, dimY); /* just copying A to U*/
        
        /* find a gradient edge-value */
        gradYsq_func(U, gradYs, dimX, dimY);
        gradval = 0.0f;
        for(i=0; i<dimY*dimX; i++) if (gradYs[i] > gradval) gradval = gradYs[i];
        gradval = gradval*gradMultip; /* max value in dy^2 gradient multipled by gradMultip */
        printf("%s %f \n", "Threshold for the sinogram is:", edgeval);
        printf("%s %f \n", "Threshold for the gradient is:", gradval);
        
        /* prepare weighting vectors */

        /* InvEucDist = (float*) calloc(pad_lookup,sizeof(float)); */ /* inverted euclidian distance (similarity) */
        /* SumInvEucDist = (float*) calloc(pad_lookup,sizeof(float)); */ /* precalculated normilization factors based InvEucDist */
        /*sumW = 0.0f;
        for(i=1; i<pad_lookup; i++) {
			InvEucDist[i] = 1.0f/sqrt(i);
			sumW += InvEucDist[i];
			SumInvEucDist[i] = sumW;
		}    
        */      
      
        /*Iterations run here*/
        iterations_number = (int)(0.25f*dimY); /*MAX iteration value cannot be larger dimY */
        /*iterations_number = 70;*/
        for(iter=0; iter < iterations_number; iter++) {
            /* calculate the gradient along Y axis to establish edges */
            if (iter > 0) {gradYsq_func(U, gradYs, dimX, dimY);}            
            /* do padding by considering all neighbours */
            pad_func_joint(U, B, gradYs, dimX, dimY, edgeval, gradval, pad_lookup, mid_val);     
            /* copy B to U*/
            copyAr(B, U, dimX, dimY);
            /* clean some outliers */
            cleaning_outliers_func(U, B, dimX, dimY, edgeval);
            /* copy B to U*/
            copyAr(B, U, dimX, dimY);
            /* check if need to terminate iterations earlier */
            counterElements = 0;
            for(i=0; i<dimX; i++) {
                for(j=0; j<dimY; j++) {
                    if (B[i*dimY + j] > edgeval) counterElements++;
                }}
            if (counterElements == dimX*dimY) {
                printf("%s \n", "Padding completed!");
                break;
            }
        }
        /* do unrolling or just copying */
        if (roll_value > 0) copyAr_roll(U, B, dimX, dimY, roll_value, 1);       
        else copyAr(U, B, dimX, dimY);         
        
        printf("%s %i \n", "Iterations stopped at:", iter);
    
    free(U);
}

float gradYsq_func(float *U, float *gradYs, int dimX, int dimY)
{
    int i,j,j1,j0;
#pragma omp parallel for shared(U, gradYs) private(i, j, j0, j1)
    for (i=0; i<dimX; i++) {
        for (j=0; j<dimY; j++) {
            j0 = j-1; if (j0 < 0) j0 = j+1;
            j1 = j+1; if (j1 == dimY) j1 = j-1;
            gradYs[i*dimY + j] = pow((0.5f*(U[i*dimY + j1] - U[i*dimY + j0])),2);
        }}
    return *gradYs;
}

float pad_func_joint(float *A, float *B, float *gradYs, int dimX,  int dimY, float edgeval, float gradval, int pad_lookup, int mid_val)
{
    int i, j, l,counter,counter1,counter2,counter3,counter4,counter5,counter6,counter7,counter8,sumcounter;
    float val1,val2,val3,val4,val5,val6,val7,val8,temp;
    
#pragma omp parallel for shared(A, B, gradYs) private(i, j, l, sumcounter, counter, val1,val2,val3,val4,val5,val6,val7,val8,counter1,counter2,counter3,counter4,counter5,counter6,counter7,counter8,temp)
    for(i=0; i<dimX; i++) {
        for(j=0; j<dimY; j++) {
                  
            /* starting a pass in a column-by-column fashion, find one pixel at the edge at a time */
            if (gradYs[i*dimY + j] >= gradval) {
				
		 /* standing on the edge pixel, we need to check if data lie in the vicinity of 
            the chosen neighbourhood of the size equal to pad_lookup */
            counter1 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && (A[i*dimY + (j + l)] >= edgeval)) {						
						counter1++; 	
                                sumcounter += counter1;
            }}	
            
            counter = 0; val1 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && (A[i*dimY + (j + l)] >= edgeval)) {
            temp = (float)(counter1 - counter)/(float)(sumcounter);
            val1 += temp * A[i*dimY + (j + l)];             
            counter++;            
            }}
            /*printf("%i %f \n", counter1, val1);*/

            counter2 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && (A[i*dimY + (j - l)] >= edgeval)) {						
						counter2++; 	
                                sumcounter += counter2;
            }}	
            
            counter = 0; val2 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && (A[i*dimY + (j - l)] >= edgeval)) {
            temp = (float)(counter2 - counter)/(float)(sumcounter);
            val2 += temp * A[i*dimY + (j - l)]; 
            counter++;            
            }}
            /*printf("%i %f \n", counter2, val2);*/

            counter3 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((i - l) >= 0) && (A[(i-l)*dimY + (j)] >= edgeval)) {						
						counter3++; 	
                                sumcounter += counter3;
            }}	
            
            counter = 0; val3 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((i - l) >= 0) && (A[(i-l)*dimY + (j)] >= edgeval)) {
            temp = (float)(counter3 - counter)/(float)(sumcounter);
            val3 += temp * A[(i-l)*dimY + (j)]; 
            counter++;            
            }}    
            /*printf("%i %f \n", counter3, val3);*/

            counter4 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((i + l) < dimX) && (A[(i+l)*dimY + (j)] >= edgeval)) {					
						counter4++; 	
                                sumcounter += counter4;
            }}	
            
            counter = 0; val4 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((i + l) < dimX) && (A[(i+l)*dimY + (j)] >= edgeval)) {
            temp = (float)(counter4 - counter)/(float)(sumcounter);
            val4 +=  temp * A[(i+l)*dimY + (j)]; 
            counter++;            
            }}
            /*printf("%i %f \n", counter4, val4);*/
   
            /* diagonal values */
            counter5 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && ((i + l) < dimX) && (A[(i + l)*dimY + (j + l)] >= edgeval)) {
						counter5++; 	
                                sumcounter += counter5;
            }}	
            
            counter = 0; val5 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && ((i + l) < dimX) && (A[(i + l)*dimY + (j + l)] >= edgeval)) {
            temp = (float)(counter5 - counter)/(float)(sumcounter);
            val5 += temp * A[(i + l)*dimY + (j + l)]; 
            counter++;            
            }}    
            /*printf("%i %f \n", counter5, val5);*/

            counter6 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && ((i - l) >= 0) && (A[(i - l)*dimY + (j + l)] >= edgeval)) {
						counter6++; 	
                                sumcounter += counter6;
            }}	
            
            counter = 0; val6 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j + l) < dimY) && ((i - l) >= 0) && (A[(i - l)*dimY + (j + l)] >= edgeval)) {
            temp = (float)(counter6 - counter)/(float)(sumcounter);
            val6 += temp * A[(i - l)*dimY + (j + l)]; 
            counter++;            
            }}  
            /*printf("%i %f \n", counter6, val6);*/

            counter7 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && ((i - l) >= 0) && (A[(i - l)*dimY + (j - l)] >= edgeval)) {
						counter7++; 	
                                sumcounter += counter7;
            }}	
            
            counter = 0; val7 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && ((i - l) >= 0) && (A[(i - l)*dimY + (j - l)] >= edgeval)) {
            temp = (float)(counter7 - counter)/(float)(sumcounter);
            val7 += temp * A[(i - l)*dimY + (j - l)]; 
            counter++;            
            }}  
            /*printf("%i %f \n", counter7, val7);*/

            counter8 = 0; sumcounter = 0;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && ((i + l) < dimX) && (A[(i + l)*dimY + (j - l)] >= edgeval)) {
						counter8++; 	
                                sumcounter += counter8;
            }}	
            
            counter = 0; val8 = 0.0f;
            for(l=1; l<pad_lookup; l++) {
    			if (((j - l) >= 0) && ((i + l) < dimX) && (A[(i + l)*dimY + (j - l)] >= edgeval)) {
            temp = (float)(counter8 - counter)/(float)(sumcounter);
            val8 += temp * A[(i + l)*dimY + (j - l)]; 
            counter++;            
            }} 		
            /*printf("%i %f \n", counter8, val8);*/

            /* total contribution of neighbors in the current position */         
            sumcounter = counter1 + counter2 + counter3 + counter4 + counter5 + counter6 + counter7 + counter8;
            /* resulting value */                        
            if (sumcounter != 0) {
            temp = (float)(counter1)/(float)(sumcounter);
            B[i*dimY + j] = temp*val1;
            temp = (float)(counter2)/(float)(sumcounter);
            B[i*dimY + j] += temp*val2;
            temp = (float)(counter3)/(float)(sumcounter);
            B[i*dimY + j] += temp*val3;
            temp = (float)(counter4)/(float)(sumcounter);
            B[i*dimY + j] += temp*val4;
            temp = (float)(counter5)/(float)(sumcounter);
            B[i*dimY + j] += temp*val5;
            temp = (float)(counter6)/(float)(sumcounter);
            B[i*dimY + j] += temp*val6;
            temp = (float)(counter7)/(float)(sumcounter);
            B[i*dimY + j] += temp*val7;
            temp = (float)(counter8)/(float)(sumcounter);
            B[i*dimY + j] += temp*val8;   }
            else  B[i*dimY + j] =  A[i*dimY + j];   
        						
        }  /* if (gradYs[i*dimY + j] >= gradval) */   
        else   B[i*dimY + j] =  A[i*dimY + j];        
       }
    }
    return *B;
}

float cleaning_outliers_func(float *A, float *B, int dimX, int dimY, float edgeval)
{
    int i, j;
#pragma omp parallel for shared(A, B) private(i, j)
    for(i=0; i<dimX; i++) {
        for(j=0; j<dimY; j++) {
            if ((j == 0) && (A[i*dimY + (j)] < edgeval)) {
                /* check the neighbour bellow */
                if (A[i*dimY + (j+1)] > edgeval)  B[i*dimY + (j)] = A[i*dimY + (j+1)];
            }
            else if ((j == (dimY-1)) && (A[i*dimY + (dimY-1)] < edgeval)) {
                /* check the neighbour above */
                if (A[i*dimY + (dimY-2)] > edgeval)  B[i*dimY + (j)] = A[i*dimY + (dimY-2)];
            }
            else {
                /* check the neighbours above and below */
                if (A[i*dimY + (j)] < edgeval) {
                    if ((A[i*dimY + (j-1)] > edgeval) && (A[i*dimY + (j+1)] > edgeval))  {
                        B[i*dimY + (j)] =  0.5*(A[i*dimY + (j-1)] + A[i*dimY + (j+1)]);   }
                }
            }
        }}
    return *B;
}

float copyAr(float *B, float *U, int dimX,  int dimY)
{
    int i;
#pragma omp parallel for shared(U, B) private(i)
    for (i=0; i<dimY*dimX; i++) {
        U[i] = B[i];
    }
    return *U;
}

float copyAr_roll(float *B, float *U, int dimX,  int dimY, int roll_value, int switcher)
{
    int i, j;
#pragma omp parallel for shared(U, B) private(i,j)
    for (i=0; i<dimX; i++) {
        for (j=0; j<dimY; j++) {
            if (switcher == 0) {
                if (j < (dimY - roll_value)) U[i*dimY + j] = B[i*dimY + (j+roll_value)];
                else U[i*dimY + j] = B[i*dimY + (j - (dimY - roll_value))];
            }
            else {
                if (j < roll_value) U[i*dimY + j] = B[i*dimY + (j+(dimY - roll_value))];
                else U[i*dimY + j] = B[i*dimY + (j - roll_value)];
            }
        }}
    return *U;
}