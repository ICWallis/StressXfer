# =============================================================
# Notes on conversion to Python:
'''
Key Changes and Notes:
Library Choices:
    - NumPy is used to handle arrays and mathematical operations.
    - Pandas is used for data handling and CSV operations.
    - matplotlib is used for plotting and saving figures.

External Functions: aphi_to_stressmags, Plot3DStressModel, and PlotStressAxes were assumed to be external functions. You'll need to implement or import these in your environment.

Table Operations: The MATLAB array2table is replaced by Pandas DataFrames, and the writetable is replaced by the to_csv method.

String Manipulations and Plotting: MATLABs scatter and text were translated using corresponding matplotlib functions.

'''
# =============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def BuildStressModel(S3_Azimuth, Aphi, Z_Obs, Z_Fault, InputElevation, MeanSurfaceElevation, path, FaultFileString):
    '''
    Build a stress model for stress transfer modeling.

    Parameters
    ----------
    S3_Azimuth : float
        The azimuth (0-360 clockwise) of the minimum horizontal stress in MPa.
    Aphi : float
        A numerical indicator of the faulting regime. 0.5 is normal faulting, 
        1.5 is strike-slip faulting, 2.5 is reverse faulting. See Simpson (1997) for more on Aphi.
    Z_Obs : array-like
        The Z values of the observation points. This value comes from the 'LoadData' function.
    Z_Fault : array-like
        The Z values of the fault facet midpoints. This value comes from the 'LoadData' function.
    InputElevation : float
        The elevation for which the initial vertical stress (Sv) is calculated, in meters.
    MeanSurfaceElevation : float
        The mean surface elevation in meters.
    path : str
        The path to the folder containing the 'RunAll.m' function.
    FaultFileString : str
        A string used to append to the file names of saved figures and exported data files.

    Returns
    -------
    Sxx_Fault, Syy_Fault, Szz_Fault, Sxy_Fault, Sxz_Fault, Syz_Fault : numpy.ndarray
        Column vectors containing the magnitudes of the calculated Cartesian stress tensor for each fault facet.
    Sxx_Obs_In, Syy_Obs_In, Szz_Obs_In, Sxy_Obs_In, Sxz_Obs_In, Syz_Obs_In : numpy.ndarray
        Column vectors containing the magnitudes of the calculated Cartesian stress tensor for each observation point.
    mu : float
        Shear modulus, relates shear stress to shear strain.
    nu : float
        Poisson's ratio.
    Myu : float
        Coefficient of friction.
    Co : float
        Cohesive strength.
    '''

    # Define constants
    ym = 75000  # Young's modulus Mpa (E)
    nu = 0.25   # Poisson's ratio
    mu = ym / (2 * (1 + nu))  # Shear modulus
    Myu = 0.60  # Coefficient of friction
    Co = 0.0    # Cohesive strength
    
    # Calculate depth in meters
    Depth = MeanSurfaceElevation - InputElevation  # Depth in meters#
    
    # Define Aphi values for plotting
    AphiAll = np.linspace(0, 3, 1001)
    
    # Call external function (assumed) 'aphi_to_stressmags'
    Shmax_MagAll, Shmin_MagAll, Sv_MagAll, _ = aphi_to_stressmags(AphiAll, Depth / 1000, Myu)
    
    L = len(AphiAll)
    Sv_MagAll = np.ones(L) * Sv_MagAll
    
    # Calculate stresses for input Aphi
    S2_Magnitude, S3_Magnitude, S1_Magnitude, _ = aphi_to_stressmags(Aphi, Depth / 1000, Myu)
    
    # Convert magnitudes to strings
    s1 = str(S1_Magnitude)
    s2 = str(S2_Magnitude)
    s3 = str(S3_Magnitude)
    
    # Pore pressure calculation
    Pp1 = (0.433 * 3280) / 145.038
    PpAll = np.ones(L) * Pp1
    
    # Plot and save stress vs Aphi figure
    NormalFaulting_InitialStressValues = plt.figure()
    plt.scatter(AphiAll, Sv_MagAll, s=2, c='r', label='Sv')
    plt.scatter(AphiAll, Shmax_MagAll, s=2, c='g', label='SHmax')
    plt.scatter(AphiAll, Shmin_MagAll, s=2, c='b', label='Shmin')
    plt.scatter(AphiAll, PpAll, s=2, c='k', label='Pp')
    plt.scatter(Aphi, S1_Magnitude, s=25, c='r')
    plt.scatter(Aphi, S2_Magnitude, s=25, c='g')
    plt.scatter(Aphi, S3_Magnitude, s=25, c='b')
    plt.scatter(Aphi, Pp1, s=25, c='k')
    
    plt.title(f'Stress at {Depth} m depth; normal faulting')
    plt.xlabel('Aphi')
    plt.ylabel('Stress (MPa)')
    plt.legend(loc='southeast')
    plt.xlim([0, 3])
    plt.ylim([0, 85])
    plt.text(0.02, 80, f'Sv = {round(S1_Magnitude, 2)}', color='red')
    plt.text(0.55, 80, f'SHmax = {round(S2_Magnitude, 2)}', color='green')
    plt.text(1.25, 80, f'Shmin = {round(S3_Magnitude, 2)}', color='blue')
    plt.text(1.9, 80, f'Pp = {round(Pp1, 2)}', color='black')
    
    # Save the figure
    filename = f'{FaultFileString}InitialStressValues.png'
    plt.savefig(f'OutputFigures/{filename}')
    
    # Build 3D stress model
    Sv_Fault, Shmax_Fault, Shmin_Fault, Sv_Obs, Shmax_Obs, Shmin_Obs, PpObs, Pp_Fault = Plot3DStressModel(
        InputElevation, Z_Obs, S1_Magnitude, S2_Magnitude, S3_Magnitude, Pp1, Z_Fault, MeanSurfaceElevation, FaultFileString
    )
    
    # Save stress model data
    StressModel = pd.DataFrame({
        'Z_Fault': Z_Fault,
        'Sv_Fault': Sv_Fault,
        'Shmax_Fault': Shmax_Fault,
        'Shmin_Fault': Shmin_Fault,
        'Pp_Fault': Pp_Fault
    })
    StressModel.to_csv(f'OutputData/{FaultFileString}StressModel.csv', index=False)
    
    # Plot stress axes for faults and observation points
    Sxx_Obs_In, Syy_Obs_In, Szz_Obs_In, Sxy_Obs_In, Sxz_Obs_In, Syz_Obs_In = PlotStressAxes(
        Sv_Obs, Shmax_Obs, Shmin_Obs, S3_Azimuth, InputElevation, FaultFileString
    )
    
    Sxx_Fault, Syy_Fault, Szz_Fault, Sxy_Fault, Sxz_Fault, Syz_Fault = PlotStressAxes(
        Sv_Fault, Shmax_Fault, Shmin_Fault, S3_Azimuth, InputElevation, FaultFileString
    )
    
    # Create trend/plunge table and save it
    S1_TP = pd.DataFrame({'Trend': [0], 'Plunge': [-90]})
    S2_TP = pd.DataFrame({'Trend': [S3_Azimuth - 90], 'Plunge': [0]})
    S3_TP = pd.DataFrame({'Trend': [S3_Azimuth], 'Plunge': [0]})
    
    # Adjust Trend if necessary
    if S3_Azimuth - 90 < 0:
        S2_TP['Trend'] = S3_Azimuth + 270
    
    # Save the trend/plunge data
    S1_TP.to_csv(f'OutputData/{FaultFileString}S1_TP_Input.csv', index=False)
    S2_TP.to_csv(f'OutputData/{FaultFileString}S2_TP_Input.csv', index=False)
    S3_TP.to_csv(f'OutputData/{FaultFileString}S3_TP_Input.csv', index=False)

# Placeholder functions for external dependencies
def aphi_to_stressmags(Aphi, Depth, Myu):
    # Implement or import aphi_to_stressmags function here
    pass

def Plot3DStressModel(InputElevation, Z_Obs, S1_Magnitude, S2_Magnitude, S3_Magnitude, Pp1, Z_Fault, MeanSurfaceElevation, FaultFileString):
    # Implement or import Plot3DStressModel function here
    pass

def PlotStressAxes(Sv, Shmax, Shmin, S3_Azimuth, InputElevation, FaultFileString):
    # Implement or import PlotStressAxes function here
    pass
