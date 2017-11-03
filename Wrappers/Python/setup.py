import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import os
import numpy
import platform	
import sys

cil_version=os.environ['CIL_VERSION']
if  cil_version == '':
    print("Please set the environmental variable CIL_VERSION")
    sys.exit(1)

extra_include_dirs = [numpy.get_include(), '../functions/']
extra_library_dirs = []
extra_compile_args = []
extra_link_args    = []
extra_libraries    = []

if platform.system() == 'Windows':
    extra_compile_args += ['/DWIN32', '/openmp']
else:
    extra_compile_args += ['-fopenmp', '-O2', '-Wall', '-std=c99']
    extra_libraries += ['m','gomp']  
    
extensions = [
    Extension("ccpi.preprocessing.limited_angle.padding", 
              sources = ["src/LimitedAngleCT/padding.pyx","../../Core/LimitedAngleCT/SinoConePad.c"],
              include_dirs = extra_include_dirs,
              library_dirs = extra_library_dirs,
              extra_compile_args = extra_compile_args,
              libraries = extra_libraries,
              extra_link_args = extra_link_args),
    ]    
setup(
     name = 'ccpi-preprocessing',
	 version = cil_version,
	 description = 'CCPi Core Imaging Library - Preprocessing Module',
	 author='Dr. Ronald Fowler',
	 license='GPL',
	 packages=['ccpi','ccpi.preprocessing','ccpi.preprocessing.beamhardening','ccpi.preprocessing.limited_angle'],
     install_requires=['numpy','scipy','matplotlib'],
	zip_safe = False,
     ext_modules = cythonize(extensions),
     package_data={'ccpi.preprocessing.beamhardening':['data/xcom/*.txt','data/spectra/Mo/*.spc','data/carouselData/*',]},
	 entry_points={
					'console_scripts': [
										'CarouselFit=ccpi.preprocessing.beamhardening.runCarouselFit:main'
									   ],
				  },
	 )
