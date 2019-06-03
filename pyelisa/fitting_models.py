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


import pandas as pd
import numpy as np
import scipy.optimize

class Lin_model(object):
    '''
    Fits data with a simple linear function.
    p1 = slope
    p2 = y intersection
    '''

    def fun(self, x, p1, p2):
        return x*p1 + p2

    def fit(self, df):
        x = []
        y = []

        for c in df:
            x.extend(df.index.values)
            y.extend(df[c].values)

        return scipy.optimize.curve_fit(self.fun, x, y, p0=(1, 0), bounds=((0, -1), (10, 3)))

    def solve(self, y, p1, p2):
        return (y-p2) / p1

class Sig_4P_model(object):
    '''
    Fits data with a sigmoidal model on a linear scale.
    
    Vector of parameters p
    p1 = Emin
    p2 = Emax
    p3 = EC50
    p4 = n_hill
    '''
    
    def fun(self, x, p1, p2, p3, p4):
        return p1 + ((p2 - p1) / (1 + (p3 / x)**p4))

    def fit(self, df):
        x = []
        y = []

        for c in df:
            x.extend(df.index.values)
            y.extend(df[c].values)

        return scipy.optimize.curve_fit(self.fun, x, y, p0=(0, 3, 100, 1), bounds=((0, 0, 0.1, 0), (1, 4, 1000, 10)))

    def solve(self, y, p1, p2, p3, p4):
        return p3 / (((p2-p1) / (y-p1)) - 1)**(1/p4)

class Sig_5P_model(object):
    '''
    Fits data with a sigmoidal model on a linear scale.
    
    Vector of parameters p
    p1 = Emin
    p2 = Emax
    p3 = EC50
    p4 = n_hill
    p5 = Asymetry factor
    '''
    
    def fun(self, x, p1, p2, p3, p4, p5):
        return p1 + ((p2 - p1) / (1 + (p3 / x)**p4)**p5)

    def fit(self, df):
        x = []
        y = []

        for c in df:
            x.extend(df.index.values)
            y.extend(df[c].values)

        return scipy.optimize.curve_fit(self.fun, x, y, p0=(0, 3, 100, 1, 1), bounds=((0, 0, 0.1, 0, 0.1), (1, 4, 1000, 10, 10)))

    def solve(self, y, p1, p2, p3, p4, p5):
        return p3 / (((p2-p1) / (y-p1))**(1/p5) - 1)**(1/p4)
