#!/bin/bash

module load R/3.6.1
srun --mpi=pmix Rscript dhmm_models.R

