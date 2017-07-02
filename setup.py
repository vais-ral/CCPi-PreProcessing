from setuptools import setup, find_packages

setup(
     name = 'ccpi-preprocessing',
	 version = '0.9',
	 description = 'CCPi Core Imaging Library - Preprocessing Module',
	 author='Dr. Ronald Fowler',
	 license='GPL',
	 packages=['ccpi','ccpi.preprocessing'],
     install_requires=['numpy','scipy','matplotlib'],
     package_dir={'':'src'},
	zip_safe = False,
     package_data={'ccpi.preprocessing':['data/xcom/*.txt','data/spectra/Mo/*.spc','data/carouselData/*',]},
	 entry_points={
					'console_scripts': [
										'CarouselFit=ccpi.preprocessing.runCarouselFit:main'
									   ],
				  },
	 )
