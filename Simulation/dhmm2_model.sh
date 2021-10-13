#!/bin/bash

module load R/3.6.1
srun --mpi=pmix Rscript dhmm2_models.R

