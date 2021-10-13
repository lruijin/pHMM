#!/bin/bash

module load R/3.6.1
srun --mpi=pmix Rscript mhmm_models.R

