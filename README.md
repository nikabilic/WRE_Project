# WRE_Project
Water Ressource Engineering Project - 2022

## PART 1: 
This week, we shall develop the aforementioned hydrological model.
Build a function (Matlab: external function; Python: use def method) that takes as inputs: our free parameters, the hourly precipitation, the monthly crop coefficient, the number of years that should be processed, and any other relevant parameter that you deem necessary (like the monthly potential evapotranspiration). 
This function should run the hydrological model for every hour, every day in a month, every month in a year and all years that should be processed.
And it should give as outputs: the generated discharge, the runoff, the infiltration, the soil saturation, the leaching, and the actual evapotranspiration.

Once you’re done, we should set up the main function that should implement the calibration (to be completed next week). This main script should:
1. read all input parameters (P, Q, etc.)
2. compute the potential evapotranspiration
3. test the hydrological model
4. compute its associated Nash-Sutcliffe coefficient by comparing it to the available discharge
