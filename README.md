# spectral_measurements

- EW.ipynb - Use EqW function from PHEW to measure equivalent widths of H alpha line for targets
- EW_measure.py - Script to run PHEW funciton on targets outside of jupyter notebook
- EqW.py - adjusted EqW function from PHEW to output EqW percentile information to .csv file
- gen_obs_tab - Use EqW function from PHEW and other tables containing terget information to generate observation table with: litname, EPIC ID, RA, DEC, local date, UTC date, UTC time, HA, sec(z), EW 16th, EW 50th, EW 84th
- overplot -function to overlay all spectal data from each exposure per target
