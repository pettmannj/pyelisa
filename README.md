# Jupyter PyELISA
Version: 1.0.0a  
Author: Johannes Pettmann
License: GNU General Public License v3

Fits ELISA standards with a 4 parameter sigmoidal function and interpolates the concentration of any unkown samples.
Open PyELISA.ipynb with Jupyter and follow the instructions.

## Requirements
* Python 3.8< (sry, I really like syntactic sugar ;-))
* ipython/jupyter
* numpy, pandas, scipy, matplotlib
* (Optional) ipywidgets for UX. See [here](https://ipywidgets.readthedocs.io/en/latest/user_install.html) for details on how to install.

### Data
Each plate data should be supplied as a CSV file. Any file in the data path will be analysed separately in one run. Files in subfolders will be ignored.

### Layout
The layout of the plate can be supplied as layout file. Layout files are CSV files that start with 'layout'. They contain the location and concentrations of the standards and blanks.
Furthermore, they can be used to adjust for any sample position. For example, if the samples are transposed when put on the plate. The results will be correctly oriented.

## Settings
Settings can be adjusted manually (e.g. setting.figure_resolution = 200) or optionally with the ipythonwidget user interface. 
If only text is shown for the cell below, ipythonwidgets are not installed correctly.
* *dilution_factor:* Dilution factor used for samples. Doesn't change the fitting. Default: 1
* *exclude_saturated:* Should oversaturated standards (i.e. OD of standards decreases at increasing doses) be excluded. Default: True
* *extrapolate_top:* How farover the top standard concentration to extrapolate. E.g. for a top standard of 500, 1.5 would allow extrapolation up to 750. Default: 1.0
* *fitting_model:* Fitting model to be used. SigModel5P (5 parameter sigmoidal model) is experimental only. Default: SigModel4P (4 parameter sigmoidal model)
* *data_path:* Folder with data files. Ignores data in subfolders. Default: /data
* *result_path:* Folder where results should be saved: Default: /results
* *data_extension:* Extension of data files. Default: .csv
* *layout_filepath:* Path to layout file. Default: layout_full.csv
* *export_results:* Export results data? Default: False
* *export_figures:* Export figures? Default: False
* *figure_type:* File format for figures. Can be .pdf, .svg, .png. Default: .pdf
* *figure_resolution:* Resolution of figure if using .png. Default: 150

**Note:** Clicking on the number next to a slider lets you change the number directly.
