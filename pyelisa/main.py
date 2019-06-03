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

import click
import yaml
import os
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as mp
import re
import csv
import pyelisa.fitting_models

def check_dir(path):
    #Check if dir exists
    if(not os.path.isdir(path)):
        print('Path "{}" does not exist'.format(path))

        #Make new folder if it doesn't exist
        try:
            os.mkdir(path)
        except:
            print('Could not create new folder')
            print('>>> Exiting application <<<')
            exit(0)

        print('New folder "{}" created'.format(os.path.split(path)[1]))

def load_settings(root_path=None, default=False):
    default_settings = {
        'model' : 3,
        'input_path' : 'data',
        'output_path' : 'results',
        'dilution_factor' : 1,
        'exclude_saturated' : True,
        'threshold_saturation' : 0,
        'extrapolate_top' : 1.1,
        'extrapolate_bottom' : 0,
        'layout' : 'layout_full.csv',
        'export_png' : False,
        'export_svg' : False
    }

    settings = None
    if(root_path == None): root_path = os.getcwd()
    path = os.path.join(root_path, 'settings.txt')
    
    if(default):
        with open(path, 'w') as f:
            yaml.dump(default_settings, f)

    try:
        with open(path, 'r') as f:
            settings =  yaml.safe_load(f)

    #If settings file is missing
    except FileNotFoundError:
        print('Settings file in "{}" not found'.format(root_path))
        print('Creating new settings file')
        load_settings(root_path, True)

        print('>>> Exiting application <<<')
        exit(0)

    #Check integrity of loaded settings (all keys from default_settings present?)
    try:
        for k in default_settings:
            settings[k]

    except KeyError:
        print('Settings file corrupted')
        print('Loading default settings')
        load_settings(True)

        print('>>> Exiting application <<<')
        exit(0)

    settings['root_path'] = root_path
    
    return settings

def load_data(settings):
    #Read data
    
    
    in_path = os.path.normpath(os.path.join(settings['root_path'], settings['input_path']))
    #out_path = os.path.normpath(os.path.join(path, settings['output_path']))

    #Load file list
    #Input files
    input_files = None
    try:
        files = os.listdir(in_path) #Get files in selected path
        input_files = [f for f in files if f[-4:] == '.csv'] #Get *.csv files in selected path
            
        #If there is no csv files
        if(not input_files):
            print('No CSV found in: "/{}"'.format(settings['input_path']))
            print('>>> Exiting application <<<')
            exit(0)

    except FileNotFoundError:
        print('Directory "/{}" does not exist'.format(settings['input_path']))
        print('Creating folder "/{}"'.format(settings['input_path']))

        try:
            os.mkdir(settings['input_path'])

        except:
            print('Could not create "/{}" directory'.format(settings['input_path']))
            print('>>> Exiting application <<<')
            exit(0)

    
    #Read CSV files
    data = {}
    
    for f in input_files:
        print(f)
        data[f[:-4]] = pd.read_csv(
            os.path.join(in_path, f),
            header=0,
            index_col=0)
    
    return data

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

def load_layout(settings):
    #Load layout file
    layout = None

    path = os.path.join(settings['root_path'], settings['layout'])

    try:
        #Read CSV file and replace empty cells (NaN) with ''
        layout = pd.read_csv(path,
            index_col=0,
            header=0,
            dtype='str').fillna('')

    except FileNotFoundError:
        print('Layout file not found')
        print('>>> Exiting application <<<')
        exit(0)

    return layout

def apply_layout(data, layout):
    blank = {}
    df_samples = {}
    df_std = {}

    pattern_std = '^[.\d]+$' #Regex for extracting standards: Start [. or Digit](1 or more) End
    pattern_samples = '^[A-H][\d][\d]?$' #Regex for extracting samples: Start [Letter A-H] [Digit][Optional 2nd digit] End
    for k in data:
        if data[k].shape != layout.shape:
            print('Error: Data file has different shape than layout file')
            print('>>> Exiting application <<<')
            exit(0)

        #Checking layout for patterns of standards and samples
        map_std = layout.applymap(regex_filter(pattern_std))
        map_samples = layout.applymap(regex_filter(pattern_samples))

        #Extract corresponding values from data and add index (well number or STD conc)
        v_std = pd.Series(data[k][map_std].values.flatten()).dropna().values
        v_samples = pd.Series(data[k][map_samples].values.flatten()).dropna().values
        idx_std = pd.Series(layout[map_std].values.flatten()).dropna()
        idx_samples = pd.Series(layout[map_samples].values.flatten()).dropna()
        std = pd.Series(v_std, index=idx_std)
        std.index = std.index.astype(float) #Convert index to float
        samples = pd.Series(v_samples, index=idx_samples)

        #Reshape sample data in correct layout (e.g. 8x10 or 8x6)
        df_samples[k] = pd.DataFrame()
        for a in samples.index:
            letter = a[0]
            number = int(a[1:])
            df_samples[k].loc[letter, number] = samples[a]

        #Reshape standards in correct layout
        df_std[k] = pd.DataFrame()
        for a in set(std.index):
            for i, v in enumerate(std[a]):
                df_std[k].loc[a, i] = v

        df_std[k].sort_index(inplace=True) #Sort standards in ascending order by index

        #Extract blank from standards and remove from standards
        blank[k] = df_std[k].loc[0].to_numpy()
        df_std[k] = df_std[k].drop(0, axis='index')

    return df_std, blank, df_samples

def export_figure(settings, title):
    '''
    Helper function for plot_fit to export figures.
    '''

    path = os.path.join(settings['root_path'], settings['output_path'])
    check_dir(path)
    
    if(settings['export_svg']):
        print('Saving figure as SVG in "{}"'.format(path))
        out = os.path.normpath(os.path.join(path, title)) #Save in export folder. Remove .csv from filename.
        mp.savefig(out + '.svg', bbox_inches='tight')

    if(settings['export_png']):
        print('Saving figure as PNG in "{}"'.format(path))
        out = os.path.normpath(os.path.join(path, title)) #Save in export folder. Remove .csv from filename.
        mp.savefig(out + '.png', bbox_inches='tight', dpi=300)

def plot_fit(df, p_fit, model, included, excluded, settings, title):
    '''
    Creates figure, showing standards, fit, included and excluded interpolated values.
    To show figure mp.show() needs to be called.
    Optionally, exports figure as PNG or SVG in root/output_path.
    
    df: Dataframe with x values as index and each y replicate in a differen column
    p_fit: Fitted parameters
    model: Instance of model used
    included: Included value as Dataframe with columns 'x' and 'y' as x and y values, respectively
    title: Title of graph
    export_svg: If SVG file should be exported
    export_png: If SVG file should be exported
    '''

    x_array = []
    y_array = []

    for c in df:
        x_array.extend(df.index.values)
        y_array.extend(df[c].values)
    
    x_fit = np.logspace(np.log10(x_array[-1]), np.log10(x_array[0]), 100)
    y_fit = model.fun(x_fit, *p_fit)

    mp.figure()
    mp.title(title)
    mp.xlabel('Standard (pg/ml)')
    mp.ylabel('OD')

    mp.semilogx(x_array, y_array, 'y.')
    mp.semilogx(x_fit, y_fit, 'k--', linewidth=0.5)
    mp.semilogx(included['x'], included['y'], '.g') #Plot included values in green
    mp.semilogx(excluded['x'], excluded['y'], '.r') #Plot excluded values in red
    
    mp.legend(['Standards', 'Fit', 'Included values', 'Excluded values'])

    #Export figure
    if(settings['export_png'] or settings['export_svg']):
        export_figure(settings, title)

def export_csv(settings, data, fit_stats):
    path = os.path.join(settings['root_path'], settings['output_path'])
    check_dir(path)
    path = os.path.join(path, 'results.csv')

    try:
        open(path, 'w').close() #Empty file
    
    except PermissionError:
        print('Permission for "results.csv" file access denied.')
        print('Make sure the file is not open in any other program.')
        print('Results will not be saved.')
        print('Please run again.')
        return
    
    write_mode = 'a' #Appending mode (adding data to the end of the file)
    
    for k in data:
        #Save Filename/dataset name
        with open(path, write_mode) as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerow([k])
            
        #Save data
        pd.DataFrame(data[k]).to_csv(path, mode=write_mode, header=None, index=None)

        #Add space and statistics
        with open(path, write_mode) as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')

            r2 = ['r2:', fit_stats[k]['r2']]
            incl = ['Included values:', fit_stats[k]['n_incl']]
            excl = ['Excluded values:', fit_stats[k]['n_excl']]
            not_fit = ['Number of not-fitted values:', fit_stats[k]['n_not_fit']]

            writer.writerow('\n')
            writer.writerow(r2)
            writer.writerow(incl)
            writer.writerow(excl)
            writer.writerow(not_fit)
            writer.writerow('\n')

    print("\n>>> Results exported <<<")

@click.command()
@click.option('-r', '--root', 'root_path', type=click.Path(True, False, True))
@click.version_option()
def main_cli(root_path):
    '''
    PyELISA  Copyright (C) 2019  Johannes
    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it under certain conditions.

    Script to automate ELISA analysis.
    It reads CSV files with OD measurments labeled with columns 1-12 and rows A-H.
    A layout files determines where to find standards and samples measurments.
    It will substract the blank, fit the standards and interpolate the samples.
    Settings for fitting can be changed, such as to exclude oversaturated standards or to define how far to extrapolate samples outside of the standard curve.

    Running the program requires passing a root path.
    Under the standard settings, all CSV files in the /data folder will be read and results exported in the /results folder.
    '''

    settings = load_settings(root_path)
    
    #Print current settings
    print('>>> SETTINGS <<<')
    for k in settings:
        print('{}: {}'.format(k, settings[k]).upper())
    print('\n')

    #Load data
    data = load_data(settings)
    layout = load_layout(settings)
    std, blank, samples = apply_layout(data, layout)
    print('Data successfully loaded\n')

    #Fitting standards
    print('Start fitting standards')
    inter_data = {} #Dict that holds the interpolated, included data
    fit_stats = {} #Holding fit stats (r2, number of included, excluded and not-fit values)
    model = None
    
    #Select fitting model
    if(settings['model'] == 1): model = fitting_models.Lin_model()
    elif(settings['model'] == 2): model = fitting_models.Sig_4P_model()
    elif(settings['model'] == 3): model = fitting_models.Sig_5P_model()
    else:
        print('Selected model not found. Please select another model.')
        print('>>> Exiting application <<<')
        exit(0)
    

    for k in data:
        #Substract blank
        std[k] = std[k] - blank[k].mean()
        samples[k] = samples[k] - blank[k].mean()

        #Exclude saturated
        if(settings['exclude_saturated']):
            max = -1
            idx = None
            for i, v in std[k].iterrows():
                #First occurance of a saturated value (i.e. not increasing anymore)
                if(v.mean() < (max * (1 + settings['threshold_saturation']))):
                    idx = std[k].loc[i:].index #Slice of this index to the end (to exclude all following values as well)
                    break

                max = v.mean()

            if(idx!=None): std[k].drop(idx, axis='index', inplace=True) #Remove saturated values
            
        #Fit standards
        popt, pcov = model.fit(std[k])

        #Calculate r2
        fit_stats[k] = {}
        y_data = []
        x_data = []

        for c in std[k]:
            x_data.extend(std[k].index.values)
            y_data.extend(std[k][c].values)

        residuals = y_data - model.fun(x_data, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data-np.mean(y_data))**2)
        fit_stats[k]['r2'] = 1 - (ss_res / ss_tot)

        #Inter-/extrapolate samples
        #The data is generated twice: Once for plotting (excluded, included) and once for exporting (inter_data)
        included, excluded = pd.DataFrame(), pd.DataFrame() #These will be series (1D array) of the included and excluded values. They won't contain any NaN values.
        inter_data[k] = pd.DataFrame(model.solve(samples[k].values, *popt)) #This is the data that will be exported. It will contain NaN values and be in the correct shape of the input (e.g. 8x10).
        mask_included = (inter_data[k] <= (settings['extrapolate_top'] * std[k].index.max())) & (inter_data[k] >= (settings['extrapolate_bottom'] * std[k].index.min()))
        included['x'] = inter_data[k][mask_included].values.flatten()
        included['y'] = samples[k].values.flatten()
        included.dropna(inplace=True) #Remove NaN values
        excluded['x'] = inter_data[k][~mask_included].values.flatten()
        excluded['y'] = samples[k].values.flatten()
        excluded.dropna(inplace=True) #Remove NaN values

        inter_data[k] = inter_data[k][mask_included] #Exclude values
        inter_data[k] = inter_data[k] * settings['dilution_factor'] #Apply dilution factor to values

        fit_stats[k]['n_incl'] = included.shape[0]
        fit_stats[k]['n_excl'] = excluded.shape[0]
        fit_stats[k]['n_not_fit'] = samples[k].size - inter_data[k].size
        
        
        plot_fit(std[k], popt, model, included, excluded, settings, k)

    mp.show()
    export_csv(settings, inter_data, fit_stats)