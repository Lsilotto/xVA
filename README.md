# Code guide 
Below we report a short guide of the Matlab code used for computing the results reported in the article "XVA Modelling: Validation, Performance and Model Risk Management" (https://link.springer.com/article/10.1007/s10479-023-05323-4). Please find further comments directly in the scripts.
We suggest to read sec. 2 of the article XVA Modelling: Validation, Performance and Model Risk Management in order to understand the logic of the codes.

The library allows to:
a) Calibrate the multicurve G2++ model in two variants: with time dependent volatility and with flat volatility. Furthermore, it is possible to calibrate G2++ parameters on ATM swaptions only or on the full swaptions cube. In particular:
  - main_calibration: G2++ parameters calibrated on ATM swaptions;
  - main_cubeCalibration: G2++ parameters calibrated on swaptions cube.

The relevant function to do this task are stored in: 
  - Calibration_Engine\HW2F_Gamma_Calibration --> G2++ model with time dependent volatility;
  - Calibration_Engine\HW2F_Calibration --> G2++ model with flat volatility;
  - Pricing_Engine\HW2F_Gamma_Pricing --> G2++ model with time dependent volatility;
  - Pricing_Engine\HW2F_Pricing --> G2++ model with flat volatility;

The calibrated parameters (both on ATM and cube) are already stored in the folder input\CalibratedParameters and can be loaded directly in the following main. 

b) Simulate the Mark-to-Future of a single instrument (Swap or Swaption) or a portfolio (of Swaps and Swaptions) and the related exposure in the presence of collateral, i.e. Variation Margin (VM) and Initial Margin (IM). Finally, the code calculate CVA, DVA, FVA and MVA for the collateralization hypotheses. In particular, the script is named main and it is divided into different sections that perform the following tasks:
  1) Import market data (rate curves, volatility and related shifts, etc.) from the Excel files in the folder input.
  2) Definition of calculation parameters (e.g. number of MC scenarios, grid granularity, etc.) and loading of calibrated G2++ parameters (please refer to point a.). In this section you can also choose whether to simulate also VM and IM. You can set the flags to 0 to save computational time.
  3) Import of instrument details from the Excel file in the folder input. This file is currently built to read IRS and Swaption records.
  4) Simulation of x and y processes using the HW2F_xy_simulation function (folder Pricing_Engine\HW2F_Gamma_Pricing).
  5) Pricing the instrument (or the portfolio) in each time step and MC scenario and computation of VM and IM. In particular, the relevant function is PortfolioExposure (folder Exposure_Engine\Simulation_Exposure) which returns the value of a portfolio of IRS and Swaptions (or the price of a 
     single instrument).
  6) Calculation of collateralized exposure, CVA, DVA, FVA and MVA.
