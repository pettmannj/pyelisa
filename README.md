# README PyELISA

PyELISA is a Python script designed to help analysing ELISA data. It imports 96 well plate data, loads a user-defined layout, fits the standards, interpolates the samples and exports the results.

## Running
Run by invoking '''PyELISA --root-path=[]'''
The root-path is the folder that contains the /data folder and the settings.txt file. Results will be saved in a /results folder.
The program will create a settings file with default settings if it is not found.
You can specifif the layout of your ELISA plate in the layout file.

## Settings
dilution_factor: Results are multiplied by this number.
exclude_saturated: If you want to automatically detect and exclude over-saturated standards.
export_png: If you want to export the standard curve as PNG...
export_svg: ..or SVG.
extrapolate_bottom: How far below the lowest standard do you want to extrapolate? 1 = lowest standard; 0 = no limit
extrapolate_top: How far above the highest standard do you want to extrapolate, e.g. 1.1 for 10% over the top standard.
input_path: Folder in root path that holds the data.
layout: Layout file in root folder.
model: Which model to fit. 1 = Linear fit; 2 = 4 parameter sigmoidal fit; 3 = 5 parameter sigmoidal fit
output_path: Folder in root path to export results.
threshold_saturation: Threshold for saturation. 0.1 means that if the next standard does not increase by at least 10%, it is consider saturated and all preceeding standards are excluded as well.

## Layout
Specify the layout in the layout.csv file. This should contain 12 columns named 1-12 and 8 rows named A-H.
Standards can be specified by number, e.g. 500 or 0 for blank.
Samples can be specified by the well number they are in, e.g. B4. The results will be reshaped based on these coordinates, i.e. a "B4" that on A1 in the input file will be on B4 in the output file.
