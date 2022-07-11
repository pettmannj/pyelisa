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

from dataclasses import dataclass
from typing import Type
from fitting_models import SigModel4P

VERSION = '1.0.0a'

#Class to store all settings
@dataclass
class Settings:
    version: str = VERSION
    dilution_factor: float = 1
    threshold_saturation: float = 1 # Needs to be larger than 1
    exclude_saturated: bool = True
    extrapolate_top: float = 1.0 # Top STD * this is the maximum concentration that is allowed for extrapolation.
    fitting_model: Type = SigModel4P

    #Export/import paths
    data_path: str = r'./data'
    data_extension: str = r'.csv' # Extension of data files to read.
    result_path: str = r'./results'
    layout_filepath: str = r'layout_full.csv'
    
    #Export data (settings are exported automatically if export_results or export_figures is turned on).
    export_results: bool = False
    export_figures: bool = False
    figure_type: str = '.pdf' # Supports pdf, svg, png and more.
    figure_resolution: int = 150 # Increasing this makes plotting take longer. Only relevant for non-vector formats (e.g. png).

    def __init__(self):
        pass
    
    def update(self, kws: dict):
        for k in kws:
            setattr(self, k, kws[k]) # Ignore key checks, since uninitialised fields do not exist yet