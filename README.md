# Dynamics of implied swaption volatility  surfaces and implications on life insurers under Solvency II

David Neuhäusler

Thesis project for MSc in Business Mathematics at [Ulm University](https://www.uni-ulm.de/)

## R code files

| Filename                   | Content                                                              |
| -------------------------- | ---------------------------------------------------------------------|
| index.Rmd                  | Empricial data analysis of swaption implied volatility surfaces      |
| 01_pca.Rmd                 | Pricipal component analysis of swaption implied volatility surfaces  |
| 02_analysis.Rmd            | Derivation of upward implied swaption volatility stress              |
| 03_hwm_calibration.Rmd     | Market-consistent calibration of HW2/G2++ model                      |
| 04_mc_simulation.Rmd       | Monte Carlo simulations for stylized insurance liability             |
| g2pp.Rmd                   | Helper functions for G2++ model analytic formulas                    |

Historic implied swaption volatility data in `swaption_volatilities.csv` is obtained from [Refinitiv Eikon](https://eikon.refinitiv.com/). The data set contains normal implied volatilities for 25 swaptions for the time period from 2/7/2013 to 9/9/2024.

## MATLAB code files

| Filename                                          | Content                                                                                                                                   |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------|
| g2pp_calibration_simulation_sensitvity_analysis.m | Market-consistent calibration of HW2/G2++ model, simulations of stylistic insurance liability and sensitivity analysis of G2++ parameters |
| LinearGaussian2Fcustom.m                          | Helper functions for G2++ model analytic formulas                                                                                         |