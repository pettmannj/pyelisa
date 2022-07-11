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

import ipywidgets as widgets
from IPython.display import display
import glob

def find_layout_files():
    layout_files = glob.glob('layout*.csv')
    return list(zip(map(lambda x: x[:-4], layout_files), layout_files))

def create_widgets():
    uxwidgets = dict()

    uxwidgets['dilution_factor'] = widgets.FloatSlider(
    value=1, # default value
    min=1, # min value
    max = 50,
    step=1, # incriment size
    description='Dilution factor:' # slider label
    )

    uxwidgets['extrapolate_top'] = widgets.FloatSlider(
    value=1, # default value
    min=1, # min value
    max = 10,
    step=0.1, # incriment size
    description='Extrapolation threshold (1 = no extrapolation):' # slider label
    )

    uxwidgets['exclude_saturated'] = widgets.Checkbox(
        value=True,
        description='Exclude saturated standards automatically?',
        indent=False
    )

    uxwidgets['export_figures'] = widgets.Checkbox(
        value=True,
        description='Export figures?',
        indent=False
    )

    uxwidgets['export_results'] = widgets.Checkbox(
        value=True,
        description='Export results (interpolated concentrations)?',
        indent=False
    )

    layout_files = find_layout_files()

    uxwidgets['layout_filepath'] = widgets.Dropdown(
        options = layout_files,
        description = 'Layout:',
        value = 'layout_full.csv' if 'layout_full.csv' in list(zip(*layout_files))[1] else layout_files[0][1] # Set layout_full as standard if present
    )

    return uxwidgets

def display_widgets(widgets=None):
    widgets = widgets if widgets else create_widgets()
    display(*widgets.values())
    return widgets

def read_widget_values(widgets):
    widget_values = dict()

    for k in widgets:
        widget_values[k] = widgets[k].value

    return widget_values