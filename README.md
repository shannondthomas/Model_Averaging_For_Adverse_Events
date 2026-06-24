# Model_Averaging_For_Adverse_Events

This repository provides all of the code associated with the paper "Identifying Treatment Effect on Adverse Events Using BMA or Stacking." The goal of these analyses was to simulation performance of model averaging methods for detecting treatment effects on adverse events across a wide range of clinical trial scenarios and provide example code using real world data. The files in this repo are as follows.

- simulation.R: Run the full simulation and output resulting data sets
- simulation_summary.R: process simulation results, compute performance metrics, compare simulated methods, run CoxPH power calculations, and generate plots/tables for publication.
- case_study_simulations.R: run simulations based on case study data to generate calibrated weight cutoff values
- case_study_runMA.R: create TTE data using the data available on Project Data Sphere, run all model averaging methods, and output results

- sim_results: contains .RData and .rds files from the main simulation
- Figures: contains all publication figures
- calibrated_results_best_table_wCoxPH.html: table of Pareto-optimal methods for all of the main simulation scenarios
- casestudy_sim_results: contains .rds files from the case study calibration simulations
