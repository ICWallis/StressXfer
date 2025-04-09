# StressXfer -- Stress Transfer Modelling

This package calculate the 3D slip that occurs on 3D input faults as a result of a 3D input stress model. The stress changes and dilatation that occur
as a result of the modeled slip at observation points that surround the input fault are also calculated. The functions then calculate the normal, shear, and 
Coulomb shear traction changes that occur fractures at the observation points. The fracture orientations used in this calculation are those that result in 
that largest normal traction reduction, shear traction increase, and Coulomb shear traction increase for all fractures oriented within 1 standard deviation
of the trend and plunge of the input faults and 1 standard deviation of the rake of the initial modeled slip.

The content of this repo is pre-release. For those wanting to use the method now, the Matlab version works. Conversion to Python is underway. 

The StressXfer_Matlab directory contains MATLAB functions that perform stress transfer modeling. Some of these functions were written or adapted from other published work.
When this is the case the citation/location of the original code is noted in the comments. The folder tree also contains the input data needed to run the
functions. Note that the input fault (.ts files) and observation points (.csv files) have geographical coordinates so that they can easily be plotted in GIS
software, these are synthetic data and not meant to represent any 'real world' location(s).  

To run these functions MATLAB software is needed (https://www.mathworks.com/products/matlab.html). These functions have only been tested with MATLAB R2020a and MATLAB R2022b.
To run, copy/unzip the folder structure exactly how it is. Navigate to the main 'StressTransferModeling' folder in MATLAB. Open the RunAll.m file. In line 13 of the 
RunAll.m code change the 'path' variable to the location of the main 'StressTransferModeling' folder on you computer. Type RunAll in the MATLAB prompt. The 
functions will generate 11 figures in the 'OutputFigures' folder, and 8 comma delimited text (.txt and .csv) files and one MATLAB workspace file in the 
'OutputData' folder.

The function is hard-coded to run with the 'AccZoneFaults.ts' file and 'Block Model_1000mcell_0melev.csv' observations points file.
There are several fault (.ts) files in the 'InputData' folder that can be used as examples of how to create your own. All of the .ts fault files and .csv
observation points file were made using Leapfrog Geothermal software (https://www.seequent.com/products-solutions/leapfrog-geothermal/). I did not attempt to
make similar files in other software.


Contributors: 
- Drew Siler (conception, method construction, MATLAB code)
- Irene Wallis (conversion of MATLAB code to Python)