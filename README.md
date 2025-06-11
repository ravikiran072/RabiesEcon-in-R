# Rabies Situational Awareness Reports

This repository contains international rabies situational awareness reports developed for public health partners. It supports automated reporting and analysis of rabies surveillance data, including data wrangling, risk classification, and geospatial visualization.

##Country Reports

- `Haiti/`: Includes cleaned REACT data, commune shapefiles, and surveillance summaries for Haiti.
- `Vietnam/`: Includes NIHE/DAH IBCM data and province-level reporting for Vietnam.



## 🔍 Script Overview for Haiti

The `script/` directory contains two core components of the monthly reporting pipeline:

- **`haiti_data_cleaning_pipeline.R`**  
  Handles data cleaning and preprocessing. It loads raw REACT surveillance data, performs spatial joins, recodes variables, and prepares epidemiological summaries for reporting.

- **`haiti_monthly_SA_report.qmd`**  
  A parameterized [Quarto](https://quarto.org/) document that generates reproducible monthly situational awareness reports. It includes dynamic maps, time series plots, and summary statistics tailored to Haiti’s IBCM surveillance data.

