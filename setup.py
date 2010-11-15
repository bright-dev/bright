#!/usr/bin/env python
 
from distutils.core import setup
#from distutils.extension import Extension

setup(name="CHAR",
	version='0.15',
	description='Cross-section Have Awesome Rates',
	author='Anthony Scopatz',
	author_email='scopatz@gmail.com',
	url='http://www.scopatz.com/',
	packages=['char'],
	package_dir={'char': 'src'}, 
	package_data={'char': ['templates/*', ' *.txt']},
    scripts=['src/scripts/char'],
	)

