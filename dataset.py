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

from typing import Dict, Type
from fitting_models import SigmoidalFittingModel
import pandas as pd
import numpy as np
import matplotlib.pyplot as mp
import settings
import os
import re
import csv
from uuid import uuid4
from datetime import date

def regex_filter(pattern):
    def fun(v):
        if v:
            m = re.search(pattern, v)
            if m:
                return True
            else:
                return False
        else:
            return False

    return fun

class Data:
    raw_data: pd.DataFrame

    # Processed data
    samples_data: pd.DataFrame
    standards_data: pd.DataFrame
    blanks_data: np.array # type: ignore

    blanked_samples: pd.DataFrame
    blanked_standards: pd.DataFrame

    # Fitting data
    # fitting_parameters
    model: SigmoidalFittingModel
    fit_stats: Dict
    interpolated_data: pd.DataFrame # Contains the interpolated concentrations of the samples. Includes also excluded values.
    inclusion_mask: pd.DataFrame # Which values are included based on the settings?

    @property
    def included_data(self):
        return self.interpolated_data[self.inclusion_mask]

    def __init__(self, model: Type) -> None:
        self.model = model()
        self.fit_stats = dict()

    def apply_layout(self, layout: pd.DataFrame):
        pattern_std = r'^[.\d]+$' # Regex for extracting standards: Start [. or Digit](1 or more) End
        pattern_samples = r'^[A-H][\d][\d]?$' # Regex for extracting samples: Start [Letter A-H] [Digit][Optional 2nd digit] End

        if self.raw_data.shape != layout.shape:
            print('ERROR: Data file has different shape than layout file')
            return

        #Checking layout for patterns of standards and samples
        map_std = layout.applymap(regex_filter(pattern_std)).values.flatten()
        map_samples = layout.applymap(regex_filter(pattern_samples)).values.flatten()

        #Extract corresponding values from data and add index (well number or STD conc)
        v_std = self.raw_data.values.flatten()[map_std]
        idx_std = layout.values.flatten()[map_std]

        v_samples = self.raw_data.values.flatten()[map_samples]
        idx_samples = layout.values.flatten()[map_samples]

        std = pd.Series(v_std, index=idx_std)
        std.index = std.index.astype(float) #Convert index to float
        samples = pd.Series(v_samples, index=idx_samples)

        #Reshape sample data in correct layout (e.g. 8x10 or 8x6)
        self.samples_data = pd.DataFrame()
        for a in samples.index:
            letter = a[0]
            number = int(a[1:])
            self.samples_data.loc[letter, number] = samples[a]

        #Reshape standards in correct layout
        self.standards_data = pd.DataFrame()
        for a in set(std.index):
            for i, v in enumerate(std[a]):
                self.standards_data.loc[a, i] = v

        self.standards_data.sort_index(inplace=True) #Sort standards in ascending order by index

        #Extract blank from standards and remove from standards
        self.blanks_data = self.standards_data.loc[0].to_numpy()
        self.standards_data = self.standards_data.drop(0, axis='index')

        # Check if there is a minimum number of different standard dilutions
        if self.standards_data.mean(axis=1).dropna().shape[0] < 4:
            raise ValueError(f'Found only {self.standards_data.mean(axis=1).dropna().shape[0]} standard dilutions. Requires at least 4 different standard dilutions.')

    def read_data(self, path: str):
        self.raw_data = pd.read_csv(
                path,
                header=0,
                index_col=0)

    def subtract_blank (self):
        avg_blank = self.blanks_data.mean()

        if np.isnan(avg_blank): # type: ignore
            raise ValueError('Blank values are missing.')

        #Substract blank
        self.blanked_samples = self.samples_data - avg_blank
        self.blanked_standards = self.standards_data - avg_blank

    # Exclude values below Emin or above Emax?
    def fit(self, setting: settings.Settings):
        #Fit standards
        self.model.fit(self.blanked_standards)

        #Calculate r2
        y_data = []
        x_data = []

        for c in self.blanked_standards:
            x_data.extend(self.blanked_standards.index.tolist())
            y_data.extend(self.blanked_standards[c].tolist())

        residuals = y_data - self.model.fun(x_data)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data-np.mean(y_data))**2)
        self.fit_stats['r2'] = 1 - (ss_res / ss_tot)

        #This is the data that will be exported. It will contain NaN values and be in the correct shape of the input (e.g. 8x10).
        self.interpolated_data = pd.DataFrame(self.model.solve(self.blanked_samples.values))
        self.interpolated_data = self.interpolated_data.fillna(0) # Fill NaNs with 0, since these are ODs below the blank.
        
        #Inter-/extrapolate samples
        #The data is generated twice: Once for plotting (excluded, included) and once for exporting (inter_data)
        max_limit = setting.extrapolate_top * self.blanked_standards.index.max()
        # min_limit = setting.extrapolate_bottom * self.blanked_standards.index.min()
        self.inclusion_mask = self.interpolated_data <= max_limit

        self.calc_fit_stats()

    def calc_fit_stats(self):
        # Calc sample size stats
        self.fit_stats['n_total'] = self.samples_data.count().sum()
        self.fit_stats['n_incl'] = self.interpolated_data[self.inclusion_mask].count().sum()
        self.fit_stats['n_excl'] = self.interpolated_data[~self.inclusion_mask].count().sum()
        
        # self.fit_stats['n_not_fit'] = self.blanked_samples.size - self.interpolated_data.size

    # Exclude saturated standards
    def exclude_saturated(self, threshold: float):
        max = -1
        idx = None
        for i, v in self.blanked_standards.iterrows():
            #First occurance of a saturated value (i.e. not increasing anymore)
            if(v.mean() < (max * threshold)):
                idx = self.blanked_standards.loc[i:].index #Slice of this index to the end (to exclude all following values as well)
                break

            max = v.mean()

        if(idx is not None):
            self.fit_stats['sat_std'] = self.blanked_standards.shape[0] - self.blanked_standards.drop(idx, axis='index').shape[0]
            self.blanked_standards.drop(idx, axis='index', inplace=True) #Remove saturated values

            # Check if there is still a minimum number of different standard dilutions
            if self.blanked_standards.mean(axis=1).dropna().shape[0] < 4:
                raise ValueError(f'Found only {self.blanked_standards.mean(axis=1).dropna().shape[0]} standard dilutions. Requires at least 4 different standard dilutions.')

            return True

        else:
            self.fit_stats['sat_std'] = 0
            return False


    def plot(self, title: str = ''):
        fig = mp.figure()
        mp.title(title)
        mp.xlabel('Standard (pg/ml)')
        mp.ylabel('OD')

        # Plot standards
        x_stds = np.array([v for v in self.blanked_standards.index.tolist() for i in range(self.blanked_standards.shape[1])])
        y_stds = self.blanked_standards.values.flatten()
        mp.semilogx(x_stds, y_stds, 'y.', label='Standards')

        # Plot fit
        x_fit = np.logspace(np.log10(x_stds.min()), np.log10(x_stds.max()), 100)  # type: ignore
        y_fit = self.model.fun(x_fit)
        mp.semilogx(x_fit, y_fit, 'k--', linewidth=0.5, label='Fit')

        # Plot sample data
        x_included = self.interpolated_data[self.inclusion_mask].values.flatten()
        x_excluded = self.interpolated_data[~self.inclusion_mask].values.flatten()
        y_data = self.blanked_samples.values.flatten()

        mp.semilogx(x_included, y_data, '.g', label='Included data') #Plot included values in green
        mp.semilogx(x_excluded, y_data, '.r', label='Excluded data') #Plot excluded values in red
  
        mp.legend()

        return fig


class Datasets:
    setting: settings.Settings
    layout: pd.DataFrame
    data: Dict[str, Data]
    iD: str

    def __init__(self, setting: settings.Settings) -> None:
        self.datasets = dict()
        self.setting = setting
        self.iD = str(uuid4())[:4] # Random ID to identify the run

        print(f'Run: {self.iD}')
        print(f'Dilution factor: {self.setting.dilution_factor}')
        print('\n')

    # Load layout file
    def load_layout(self):
        path = self.setting.layout_filepath

        try:
            #Read CSV file and replace empty cells (NaN) with ''
            self.layout = pd.read_csv(path,
                index_col=0,
                header=0,
                dtype='str').fillna('')

        except FileNotFoundError:
            print('ERROR: Layout file not found')

        print(f'Loaded layout file {path}.')

    def load_data(self):
        #Read data
        in_path = self.setting.data_path

        if not os.path.isdir(in_path):
            print(f'ERROR: Could not find data path: {in_path}')
            return

        #Load file list
        #Input files
        input_files = None

        files = os.listdir(in_path) #Get files in selected path
        input_files = [f for f in files if f[-4:] == self.setting.data_extension] #Get *.csv files in selected path
            
        #If there is no data files files
        if(not input_files):
            print(f'No {self.setting.data_extension} found in: "/{self.setting.data_path}"')

        #Read CSV files and apply layout
        self.data = dict()

        for f in input_files:
            print(f'Loaded data file {f}.')

            key = f[:-4]
            self.data[key] = Data(self.setting.fitting_model)
            path = os.path.join(in_path, f)

            try:
                self.data[key].read_data(path)
                self.data[key].apply_layout(self.layout)
                self.data[key].subtract_blank()

                # Exclude saturated
                if self.setting.exclude_saturated:
                    if self.data[key].exclude_saturated(self.setting.threshold_saturation):
                        print(f'Removed saturated standards for dataset {key}.')
                    
                    else:
                        print(f'No saturated standards identified for dataset {key}.')

            except ValueError as e:
                print(f'ERROR in file {f}. Please check if {self.setting.layout_filepath} is the correct layout file and that the data file is correctly formatted.')
                raise e


    def load(self):
        self.load_layout()
        self.load_data()

    def export_csv(self):
        path = os.path.join(os.path.normpath(self.setting.result_path), f'{self.iD} - results.csv')

        try:
            open(path, 'w').close() #Empty file
        
        except PermissionError as e:
            print('Permission for "results.csv" file access denied.')
            print('Make sure the file is not open in any other program.')
            print('Results will not be saved.')
            print('Please run again.')
            raise e
        
        write_mode = 'a' #Appending mode (adding data to the end of the file)
        with open(path, write_mode) as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')

            # Write header
            writer.writerow(['Run ID', self.iD])
            writer.writerow(['Date', date.today().strftime("%B %d, %Y")])
            writer.writerow(['PyELISA version', self.setting.version])
            writer.writerow(['Dilution factor', self.setting.dilution_factor])
            writer.writerow('\n\n')
            
            for k in self.data:
                d = self.data[k].included_data.copy() * self.setting.dilution_factor
                d = d.fillna('SATURATED') # Fill saturated values with string
                d['sep'] = 'X' # Add Xs on right border

                width = d.shape[1] # Number of columns
                writer.writerow([k] + (['X'] * (width-1))) # Write header with name + Xs
                d.to_csv(csvfile, header=None, index=None, line_terminator='\n')  # type: ignore
                writer.writerow(['X'] * width) # Add Xs below

                # Add stats
                stats = self.data[k].fit_stats
                writer.writerow(['Number of oversaturated standard dilutions excluded', stats.get('sat_std', 'Turned off')])
                writer.writerow(['r2', stats['r2']])
                writer.writerow(['Total number of samples', stats['n_total']])
                writer.writerow(['Interpolated samples', stats['n_incl']])
                writer.writerow(['Samples over extrapolation threshold', stats['n_excl']])
                writer.writerow('\n\n')

        print("\n>>> Results exported <<<")

    def fit(self):
        #Fitting standards
        print('\n>>> Start fitting standards <<<')
        # inter_data = {} #Dict that holds the interpolated, included data
        # fit_stats = {} #Holding fit stats (r2, number of included, excluded and not-fit values)

        for k in self.data:
            # Fit data
            self.data[k].fit(self.setting)
               
            # Plot fit
            fig = self.data[k].plot(k)

            # Save figure
            if self.setting.export_figures:
                path = os.path.normpath(os.path.join(self.setting.result_path, f'{self.iD} - {k}'))
                extension = self.setting.figure_type
                fig.savefig(path + extension, bbox_inches='tight', dpi=self.setting.figure_resolution)

        mp.show()

        if self.setting.export_results: self.export_csv()