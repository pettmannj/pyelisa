# -*- coding: utf-8 -*-

'''
PyELISA
Python script to analyse ELISA data
Copyright (C) 2019  Johannes Pettmann

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from setuptools import setup, find_packages

setup(
    name='PyELISA',
    version='0.1',
    description='A program for automated analysis of ELISA data.',
    author='Johannes Pettmann',
    license='GNU GPLv3',
    python_requires='>=3.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click >= 7.0',
        'Pandas >= 0.24',
        'NumPy >= 1.16',
        'pyyaml >= 5.1',
        'matplotlib >= 3.0',
        'scipy >= 1.2'
    ],
    entry_points = {
        'console_scripts': ['PyELISA=pyelisa.main:main_cli'],
      },
)