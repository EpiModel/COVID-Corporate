# Modeling the impact of adherence to U.S. isolation and masking guidance on SARS-CoV-2 transmission in office workplaces in 2021-2022

This repository holds the source code to reproduce the analysis of our model of SARS-CoV-2 transmission dynamics among office workers and their communities in Georgia, USA in 2021/2022.

These models are written and executed in the R statistical software language. To run these files, it is necessary to first install our epidemic modeling software, [EpiModel](https://github.com/EpiModel/EpiModel), and our extension package specifically for modeling SARS-CoV-2 transmission dynamics, [EpiModelCOVID](https://github.com/EpiModel/EpiModelCOVID). The branch of the EpiModelCOVID repository associated with this research project is [Corporate-2](https://github.com/EpiModel/EpiModelCOVID/tree/Corporate-2).

The `R` folder contains all code required to run the analysis, including code to analyze the Corporate Mix contact study to parametrize social contacts in the transmission model, and the transmission model code.
