#!/bin/bash

export OMP_NUM_THREADS=1

mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum,eig}}

# Lattice dims
LX=8
LY=8

# The value of the coupling in the U(1) 2D theory
BETA=3.0

# The total number of thermalised HMC iterations to perform.
HMC_ITER=1000
# The number of HMC iterations for thermalisation (accept + accept/reject).
HMC_THERM=10
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=5
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=50
# Reverse the gauge fields
HMC_REVERSE=100
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0

# HMC time steps in the integration 
HMC_NSTEP=3
# HMC inner time steps in the integration 
HMC_INNER_NSTEP=1
# Degree of polynomial for AlgRemez
HMC_AR_DEGREE=12
# Precision for AlgRemez (GMP)
HMC_AR_GMP_PREC=40

# HMC trajectory time
HMC_TAU=1.0
# Integrator type: leapfrog = 0, fgi = 1
# FYI, aim for 70% acceptance with Leapfrog
# and 90% with FGI for optimal FLOP usage
INTEGRATOR=1

# Number of APE smearing hits to perform when measuring topology
APE_ITER=1
# The alpha value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=1

# Dynamic fermion parameters
FLAVOURS=2
# Light Fermions (degenerate) mass
MASS=0.1
# Heavy Fermion mass
MASS_HEAVY=0.6


# Maximum CG iterations
MAX_CG_ITER=10000
# CG eps (tolerance)^2
CG_EPS=1e-18

# Eigensolver parameters
INSPECT_SPECTRUM=0
DEFLATE=0
NKR=56
NEV=24
NCONV=24

# Tolerance on the residual
EIG_TOL=1e-10
# Maximum restart iterations
MAXITER=1000

#polyACC (experimental)
USE_ACC=0
AMAX=11
AMIN=1.0
N_POLY=100

X_BLK=4
Y_BLK=4
N_LOW=16
NDEFL=${N_LOW}

# Measuremets: 1 = measure, 0 = no measure
# Polyakov loops
MEAS_PL=1
# Wilson loops and Creutz ratios
MEAS_WL=1
# Pion Correlation function
MEAS_PC=1
# Vacuum trace
MEAS_VT=1

command="./wilson2D $BETA $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT 
         $HMC_CHKPT_START $HMC_NSTEP $HMC_TAU $APE_ITER $APE_ALPHA $RNG_SEED 
	 $DYN_QUENCH $MASS $MAX_CG_ITER $CG_EPS $DEFLATE $NKR $NEV $NCONV 
	 $EIG_TOL $MAXITER $USE_ACC $AMAX $AMIN $N_POLY $X_BLK $Y_BLK $N_LOW $NDEFL 
	 $MEAS_PC $MEAS_WL $LX $LY $HMC_INNER_NSTEP $MASS_HEAVY $INTEGRATOR $FLAVOURS
	 $HMC_REVERSE $HMC_AR_DEGREE $HMC_AR_GMP_PREC $INSPECT_SPECTRUM"

echo $command

$command
