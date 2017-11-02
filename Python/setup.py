from setuptools import setup, find_packages
import os

cil_version=os.environ['CIL_VERSION']
if  cil_version == '':
    print("Please set the environmental variable CIL_VERSION")
    sys.exit(1)

setup(
     name = 'ccpi-preprocessing',
	 version = cil_version,
	 description = 'CCPi Core Imaging Library - Preprocessing Module',
	 author='Dr. Ronald Fowler',
	 license='GPL',
	 packages=['ccpi','ccpi.preprocessing','ccpi.preprocessing.beamhardening'],
     install_requires=['numpy','scipy','matplotlib'],
	zip_safe = False,
     package_data={'ccpi.preprocessing.beamhardening':['data/xcom/*.txt','data/spectra/Mo/*.spc','data/carouselData/*',]},
	 entry_points={
					'console_scripts': [
										'CarouselFit=ccpi.preprocessing.beamhardening.runCarouselFit:main'
									   ],
				  },
	 )
