# -*- coding: utf-8 -*-

'''
PyELISA
Python script to analyse ELISA data
Copyright (C) 2021  Johannes Pettmann

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
from typing import Tuple, Sequence

class SigmoidalFittingModel:
    _bounds: Tuple[Tuple[float, ...], Tuple[float, ...]]
    _initial_guess: Tuple[float, ...]
    _fitting_results: Sequence
    _emin: float = 0 # Fix Emin to this value

    def fun(self, x, *parameters):
        raise NotImplementedError

    def solve(self, y, *parameters):
        raise NotImplementedError

    def prepare_fit(self, data: pd.DataFrame):
        x = []
        y = []

        for c in data:
            x.extend(data.index.tolist())
            y.extend(data[c].tolist())

        return x, y

    def fit(self, data: pd.DataFrame):
        x, y = self.prepare_fit(data)
        self._fitting_results, _ = scipy.optimize.curve_fit(self.fun, x, y, p0=self._initial_guess, bounds=self._bounds)

    def get_emin(self):
        raise NotImplementedError

class SigModel4P(SigmoidalFittingModel):
    '''
    Fits data with a sigmoidal model on a linear scale.
    
    Vector of parameters p
    Emin (fixed to 0)
    p1 = Emax
    p2 = EC50
    p3 = n_hill
    '''

    def __init__(self) -> None:
        self._bounds = (0, 0.1, 0), (10, 1000, 10)
        self._initial_guess = 3, 100, 1
    
    def fun(self, x, *parameters):
        parameters = parameters if parameters else self._fitting_results
        return self._emin + ((parameters[0] - self._emin) / (1 + (parameters[1] / x)**parameters[2]))
        # return p1 + ((p2 - p1) / (1 + (p3 / x)**p4))

    def solve(self, y, *parameters):
        parameters = parameters if parameters else self._fitting_results

        data = np.where(y > 0, y, np.NaN) # Exclude negative data
        data = np.where(data < parameters[0], data, np.NaN) # Exclude data over the Emax
        data = np.where(data > self._emin, data, np.NaN) # type: ignore Exclude data below the Emin

        return parameters[1] / (((parameters[0] - self._emin) / (data - self._emin)) - 1)**(1/parameters[2])  # type: ignore


class SigModel5P(SigmoidalFittingModel):
    '''
    EXPERIMENTAL!!!
    Fits data with a sigmoidal model on a linear scale.
    
    Vector of parameters p
    Emin (fixed to 0)
    p1 = Emax
    p2 = EC50
    p3 = n_hill
    p4 = Asymetry factor
    '''

    def __init__(self) -> None:
        self._bounds = (0, 0, 0.1, 0, 0.1), (1, 10, 1000, 10, 10)
        self._initial_guess = 0, 3, 100, 1, 1
    
    def fun(self, x, parameters=None):
        parameters = parameters if parameters else self._fitting_results
        return parameters[0] + ((parameters[1] - parameters[0]) / (1 + (parameters[2] / x)**parameters[3])**parameters[4])

    def solve(self, y, *parameters):
        parameters = parameters if parameters else self._fitting_results

        data = np.where(y > 0, y, np.NaN) # Exclude negative data
        data = np.where(data < parameters[0], data, np.NaN) # Exclude data over the Emax
        data = np.where(data > self._emin, data, np.NaN) # type: ignore Exclude data below the Emin

        return parameters[1] / (((parameters[0] - self._emin) / (data - self._emin))**(1/parameters[3]) - 1)**(1/parameters[2])  # type: ignore
