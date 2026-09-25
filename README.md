# Windfarm-Bivalve-Habitat-Overlap

# Code for Sun et al., Disproportionate Overlap Between Offshore Wind Development Footprint and Suitable Bivalve Habitat in the Mid-Atlantic (manuscript under review).

Overview
This repository contains R scripts used to model habitat suitability for Atlantic surfclam, ocean quahog, and Atlantic sea scallop in the U.S. Mid-Atlantic and compare suitable habitat with offshore wind energy areas (WEAs). The analyses include:
Processing survey observations and environmental layers (bottom temperature, depth, and sediment grain size)
Fitting variable-specific suitability curves and combining them into habitat suitability indices using boosted regression tree weights
Evaluating habitat suitability with leave-one-year-out cross-validation
 Mapping habitat suitability and quantifying its overlap with WEAs
Classifying environmental regimes with a self-organizing map and comparing their representation inside and outside WEAs
Producing figures and summary statistics for the manuscript

The analyses use NEFSC shellfish and VIMS scallop survey data, GLORYS bottom temperature, GEBCO bathymetry, USGS sediment data, and BOEM wind energy area boundaries. Raw data and generated results are not included here; the scripts refer to local ../data and ../results directories. Obtain input data from the original providers before running the workflow.



# Citation
If you use this code, please cite the associated manuscript once a publication reference is available. For questions about the data or workflow, contact Ming Sun at ming.sun@vims.edu.
