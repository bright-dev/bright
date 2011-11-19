#!/usr/bin/env python
 
from distutils.core import setup

setup(name="CHAR",
	version='1.0',
	description='Cross-section Have Awesome Rates',
	author='Anthony Scopatz',
	author_email='scopatz@gmail.com',
	url='http://www.scopatz.com/',
	packages=['char', 'char.templates', 'char.templates.lwr', 
              'char.ui', 'char.run'],
	package_dir={'char': 'char'}, 
    scripts=['char/scripts/char'],
	)

