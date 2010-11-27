#!/usr/bin/env python
 
from distutils.core import setup
from setuptools.command.develop import develop
#from distutils.extension import Extension

setup(name="CHAR",
	version='0.15',
	description='Cross-section Have Awesome Rates',
	author='Anthony Scopatz',
	author_email='scopatz@gmail.com',
	url='http://www.scopatz.com/',
	packages=['char', 'char.templates', 'char.templates.lwr'],
	package_dir={'char': 'char'}, 
    scripts=['char/scripts/char'],
	)

