LiDAR Reveals Complex Trade-Offs Between Coppicing, Forest Structure, and Biodiversity in UK Woodland Restoration

This repository contains the data processing, analysis, and figure-generation code supporting the study:
“LiDAR Reveals Complex Trade-Offs Between Coppicing, Forest Structure, and Biodiversity in UK Woodland Restoration.”

The work uses airborne LiDAR–derived structural metrics, combined with ecological survey data, to quantify how coppicing management influences forest structure and biodiversity responses in UK woodlands.

├── code/               # Analysis scripts (see below)
├── data/               # Raw and processed data (Some sensitive data not included; see Data availability)
├── figures/            # Figures
├── manuscript/         # Manuscript files (drafts / submitted version)
├── supp_info/          # Supplementary information
├── README.md
└── .gitignore

Code overview

The code/ directory contains scripts used to:

Process airborne LiDAR data and derive forest structural metrics

Integrate LiDAR metrics with biodiversity and management data

Fit statistical models linking coppicing, structure, and biodiversity responses

Generate figures and tables used in the manuscript


Data availability

Some raw LiDAR point clouds and some ecological datasets are not included in this repository due to data size, licensing, and access restrictions.

data/ contains placeholder files and README notes describing how raw data were obtained and preprocessed. It also contains derived metrics and cleaned datasets sufficient to reproduce analyses and figures, where redistribution is permitted.

Please contact the author for guidance on accessing restricted datasets.

Reproducibility

Analyses were conducted primarily in R.
For full reproducibility, users should ensure package versions are compatible with those current at the time of analysis. An renv environment may be added in future updates.

How to run the analysis:

Clone the repository
Ensure required R packages are installed
Populate data_raw/ according to instructions provided
Run scripts in the code/ directory
Outputs (figures and tables) will be written to figures/ and outputs/

Citation

If you use this code or data, please cite the associated manuscript:

Beese, L. V. J., & Jucker, T. (2026). LiDAR Reveals Complex Trade-Offs Between Coppicing, Forest Structure, and Biodiversity in UK Woodland Restoration. Manuscript in preparation.

Author

Lucy Beese
University of Bristol

Licence

This repository is shared for academic and non-commercial research purposes.
A licence file will be added prior to final publication.
