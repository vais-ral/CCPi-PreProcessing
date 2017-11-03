
"""
Copyright 2017 CCPi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

SinoPad.pyx
A Sinogram padding for limited angle reconstruction
Author: Dr. Daniil Kazantsev 
"""
import cython

# import numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern void SinoPad(float *A, float sinoMultip, float gradMultip, int pad_lookup, int roll_value, int dimX, int dimY, float *B, float *gradYs) 

@cython.boundscheck(False)
@cython.wraparound(False)
def sinogram_pad(np.ndarray[np.float32_t, ndim=2, mode="c"] sinogram, float sinoMultip, float gradMultip, int pad_lookup, int roll_value):
	"""
	Takes an input sinogram and using the padding parameters this algorithm pads the missing angles
	and returns a sinogram.
	
	param A: singoram [detectors, angles] 
	param sinoMultip: a multiplier for the sinogram in the range (0,1];   
	param gradMultip: a multiplier for the gradient in the range (0,1];   
	param pad_lookup:  number of pixels to add to calculate a median of usabale values 
	param roll_value: symmetric sinogram roll from the top to bottom, 0 - do not roll 
	"""
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] padded_sinogram = np.empty([sinogram.shape[0], sinogram.shape[1]], dtype='float32')
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] gradients = np.empty([sinogram.shape[0], sinogram.shape[1]], dtype='float32')
	SinoPad(&sinogram[0,0], sinoMultip, gradMultip, pad_lookup, roll_value, sinogram.shape[0], sinogram.shape[1],&padded_sinogram[0,0], &gradients[0,0])
	return padded_sinogram, gradients
