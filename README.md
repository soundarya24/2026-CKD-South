This repo builds the CKD-CC-temperature association.

How to run the analysis:
1. first run the r file in data management folder, that preps the data for analysis; check whether derived data is built
2. then proceed to analysis and run the r file for association; this will populate the docs folder

This is the directory structure

.
├── 1_dataManagement
│   └── 00_build_panel.R
├── 2026-CKD-South.Rproj
├── 2_derivedData
├── 3_analysis
│   └── 01_aim1_association.R
├── 4_docs
├── README
└── _rawData
    ├── ckd-all-measures
    │   ├── citation.txt
    │   ├── daly_yll_yld.csv
    │   └── deaths_inc_prev.csv
    └── temperature
        ├── andhrapradesh_percentiles.csv
        ├── karnataka_percentiles.csv
        ├── kerala_percentiles.csv
        ├── tamilnadu_percentiles.csv
        └── telangana_percentiles.csv
