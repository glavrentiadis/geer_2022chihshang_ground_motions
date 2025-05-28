# Extracting and Characterizing Near-fault Effects from Strong Motion Observations of the 2022 Chishang Earthquake Sequence: Methodology and Validation

## Description

This repository contains the signal processing, analysis, and visualization code associated with the manuscript:

Lavrentiadis G, Asimaki D, Carey TJ, and Mason HB  *Extracting and Characterizing Near-Fault Effects from Strong Motion Observations of the 2022 Chihshang Earthquake Sequence: Methodology and Validation*.  (Currently under review at the *Bulletin of the Seismological Society of America*)


## Repository Structure

The main folder, `Analyses`, contains all the scripts for preprocessing, baseline correction, regression, validation, and supporting libraries.

Within the `Analyses` folder:
- `ground_motions/` contains scripts and subfoders for:  
  - formatting raw time histories
  - `tilt_correction/` applying tilt-based baseline correction
  - `vel_pulse/` identifying and characterizing velocity pulses
  - `summary/` summarizing ground motion intensity measures
- `gis/` contains QGIS project files and scripts for generating geospatial layers.
- `python_lib/` and `matlab_lib/` contain custom functions for evaluating velocity models.

The `Data` folder mirrors the structure of the `Analyses` folder and includes all the corresponding input and output files.

The `Raw_files` folder unprocessed project files in their original format.



```
    .
    |--Analyses
    |     |--ground_motions
    |     |    |--tilt_correction
    |     |    |--vel_pulse
    |     |    |--summary
    |     |
    |     |--gis
    |     |--python_lib
    |     |--matlab_lib
    |
    |--Data
    |     
    |--Raw_files
```

## Collaborators
 - Grigorios Lavrentiadis -- Postdoctoral Associate, California Institute of Technology
 - Domniki Asimaki -- Professor, California Institute of Technology
 - Trevor J Carey -- Assistant Professor, The University of Toronto
 - Henry Benjamin Mason -- Associate Professor, University of Nevada, Reno

## Acknowledgments 
The authors gratefully acknowledge the support of the Geotechnical Extreme Events Reconnaissance (GEER) Association. Financial support for GEER is provided by the National Science Foundation (NSF) through the Geotechnical Engineering Program under Grant No. CMMI-1266418. 
Any opinions, findings, conclusions, or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation.

