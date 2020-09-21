Everything you need to know to use this grid is in this file, including the updates.
This file lists the file structure and nomenclature of the grid.

For more details please Refer, Goyal et al. 2020,MNRAS, "A Library of Self-consistent Simulated Exoplanet Atmospheres" 

***Updates*** 

20th March 2020. 
1) Self-Consistent Grid Uploaded Online

************************************************************

Planetary Parameters are adopted from TEPCAT database (Southworth 2011a).
All these files have been compressed to reduce their size, therefore decompress them using “gzip -d filename” command to obtain .txt file.

1) Pressure-Temperature Profile files

Column 1 - Pressure (bar)
Column 2 - Temperature (K)

File Name Nomenclature: pt-eqpt_PlanetName_Recirculationfactor_log(Metallicity)_C/ORatio_model.txt.gz

A quick python plotting program for transmission spectra is also provided “quick_pt_plotting.py” in “Utilities” folder

2) Chemical Abundances Files

Column 1 - Pressure (bar)
Column 2-278 - Molecular Abundances (Mole Fraction)

File Name Nomenclature: chem-eqpt_PlanetName_Temperature_log(Metallicity)_C/ORatio_model.txt.gz

A quick python plotting program for chemical abundances is also provided “quick_chemical_abundances_plotting.py” in “Utilities” folder

(See chemistry_files_column_names.txt in “Utilities” folder for ordered list (column number) of all molecules)


3) Transmission Spectra File

Column 1 - Wavelength in microns
Column 2 - Transit Depth ((Rp/Rs)^2)

File Name Nomenclature = trans-eqpt_PlanetName_Recirculationfactor_log(Metallicity)_C/ORatio_model.txt.gz

A quick python plotting program for transmission spectra is also provided “quick_trans_spectra_plotting.py” in “Utilities” folder

4) Emission Spectra File

Header Line 1 - Planetary Radius (Rp_TOA) and Stellar Radius (Rs) are given in centimeters(cm).
Column 1 - Wavelength (microns)
Column 2 - Stellar Flux (Fs) (Watts per meter square per micrometer)
Column 3 - Planetary Flux (Fp) (Watts per meter square per micrometer)

File Name Nomenclature = emiss-eqpt_PlanetName_Recirculationfactor_log(Metallicity)_C/ORatio_model.txt.gz

A quick python plotting program for emission spectra is also provided “quick_emiss_spectra_plotting.py” in “Utilities” folder


5) Parameter files (If needed) 

These files provided in “Utilities/each_planet_grid_parameter_table/“ folder have arrays of grid parameters for each planet which can be used for copying or directly importing to Python

Parameter File Grid Points
Row 1 - Temperature (K)
Row 2 - Metallicity (logarithmic)
Row 3 - C/O ratio
Row 4 - Haze
Row 5 - Cloud


Note - Please cite Goyal et al. 2020, MNRAS, "A Library of Self-consistent Simulated Exoplanet Atmospheres"
while using this data

Contact - Jayesh Goyal (jgoyal@astro.cornell.edu)
Please contact above email address for any bugs, suggestions and clarification of this exoplanet atmospheres database.
 
 

