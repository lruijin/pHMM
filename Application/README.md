## Working flow for reproducing real data analysis results
### Obtain posterior samples
The three files in this folder `phmm_sel.R`, `pHMMe_sel.R` and `mhmm_pc_sel.R` are used to analyze the FMOD trail data with perception-augmented HMM, extended perception-augmented HMM and multivariate mixed effect HMM, respectively. Each of them consists of three steps:
 - Step 1: read in data, including manifesto variables, covarites and treatment arms.
 - Step 2: set initial values. For different number of hidden states, there are different numbers of parameters to be estimated. Thus, we different initial values are specified for 2-, 3-, 4- and 5-class models. The 3-class model is the primary model we used in the paper and the default setting in the code.
 - Step 3: obtain the posterior samples of the parameters. The final sample includes 40000 samples, with the first 30000 samples as burn-in.
### Summarize the estimation results from the posterior samples
### Model diagnosis
