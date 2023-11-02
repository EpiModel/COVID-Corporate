#!/bin/bash

. /projects/epimodel/spack/share/spack/setup-env.sh
spack load r@4.3.0

Rscript hpc_calib_vax.R
