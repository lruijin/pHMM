# Perception-augmented Hidden Markov Model
This repository contains the codes for reproducing results in "A Perception-augmented Hidden Markov Model for Family Management of Diabetes". It consists of three main components:
1. Folder `Application` contains codes for analyzing FMOD trail data and reproducing results in Section 2 and 5
2. Folder `Simulation` contains codes for reproducing results of simulation studies in Section 6
3. `DHMM_*.R` are the main functions implementing the proposed method and `MHMM_*.R` implement the method introduced in Raffa and Dubin (2015) for the purpose of comparison:
  - `DHMM_v1.R` and `MHMM_v3.R` are used for simulation study. 
  - `DHMM_v2.R` and `MHMM_v4.R` are used for real data application. They accommodate covariates in the emission model.  

The working flows for generating real data analysis results (as in Sections 2 and 5) and simulation results (as in Section 6) are described in separate `README` files in folders `Application` and `Simulation`, respectively. 
## Reference:
1. Lu, R., Nansel, T. R. and Chen, Z. (2021). "A perceptiion-augmented hidden Markov model for family management of diabetes". Under review.
2. Raffa, J. D. and Dubin, J. A. (2015). "Multivariate longitudinal data analysis with mixed effects hidden Markov models." Biometrics, 71, 3, 821-831.
