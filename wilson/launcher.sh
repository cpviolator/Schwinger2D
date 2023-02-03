#!/bin/bash

export OMP_NUM_THREADS=1

mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum,eig}}

# The RNG seed
SEED=1234

# General verbosity
VERBOSITY=true

# Physics params
#-----------------------------------------------
# Lattice dims
LX=10
LY=10
# The value of the coupling in the U(1) 2D theory
BETA=5.0
# Dynamic fermion parameters
# 0 = pure gauge
# 2 = two light degenerate fermions
# 3 = two light degenerate, one heavy fermion
FLAVOURS=2
# Light Fermions (degenerate) mass
MASS=0.01
# Heavy Fermion mass
MASS_HEAVY=0.6

# HMC params
#-----------------------------------------------
# The total number of thermalised HMC iterations to perform.
HMC_ITER=1000
# The number of HMC iterations for thermalisation (accept + accept/reject).
HMC_THERM=00
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=10
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=100
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=100
# Reverse the gauge fields for ergodicity check
HMC_REVERSE=100
# HMC time steps in the integration 
HMC_NSTEP=4
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
HMC_INTEGRATOR=1
# Maximum CG iterations
MAX_CG_ITER=10000
# CG residual tolerance 
CG_TOL=1e-9
# CG verbosity
CG_VERBOSITY=false

# Eigensolver parameters
#-----------------------------------------------
NKR=200
NEV=200
NCONV=200
NDEFL=100
# IRAM MG projector
X_BLK=2
Y_BLK=2
N_LOW=50
# Tolerance on eigenvector residua
EIG_TOL=1e-10
# Maximum restart iterations
MAX_IRAM_ITER=1000
# Operator to solve
EIG_OP=MdagM
# Spectrum to compute
EIG_SPEC=SR
# IRAM verbosity
EIG_VERBOSITY=false
# inspect spectrum
EIG_INSPECTION=true

# Measurements
#-----------------------------------------------
# Number of APE smearing hits to perform when measuring topology
APE_ITER=0
# The alpha value in the APE smearing
APE_ALPHA=0.5
# Measurements: 1 = measure, 0 = no measure
# Polyakov loops
MEAS_PL=0
# Wilson loops and Creutz ratios
MEAS_WL=0
# Pion Correlation function
MEAS_PC=0
# Vacuum trace
MEAS_VT=0

# Execute
#-----------------------------------------------


BASIC_PARAMS="--seed ${SEED} --verbosity ${VERBOSITY} --beta ${BETA} --dim ${LX} ${LY} \
	      --mass ${MASS} --mass-heavy ${MASS_HEAVY} --flavours ${FLAVOURS} "

HMC_PARAMS="--hmc-traj-length ${HMC_TAU} --hmc-n-step ${HMC_NSTEP} --hmc-inner-step ${HMC_INNER_NSTEP} \
	    --hmc-n-trajectories ${HMC_ITER} --hmc-therm ${HMC_THERM} --hmc-checkpoint ${HMC_CHKPT} \
	    --hmc-checkpoint-start ${HMC_CHKPT_START} --hmc-skip ${HMC_SKIP} --hmc-reverse ${HMC_REVERSE} \
            --hmc-integrator ${HMC_INTEGRATOR} \
	    --pfe-degree ${HMC_AR_DEGREE} --pfe-prec ${HMC_AR_GMP_PREC} \
	    --cg-max-iter ${MAX_CG_ITER} --cg-tol ${CG_TOL} --cg-verbosity ${CG_VERBOSITY}"

EIG_PARAMS="--eig-n-ev ${NEV} --eig-n-kr ${NKR} --eig-n-conv ${NCONV} --eig-n-deflate ${NDEFL} --eig-max-restarts ${MAX_IRAM_ITER} \
            --eig-tol ${EIG_TOL} --eig-operator ${EIG_OP} --eig-spectrum ${EIG_SPEC} --eig-block-scheme ${X_BLK} ${Y_BLK} \
            --eig-low-modes ${N_LOW} --eig-verbosity ${EIG_VERBOSITY} --eig-inspection ${EIG_INSPECTION}"

MEASUREMENTS="--ape-alpha ${APE_ALPHA} --ape-iter ${APE_ITER} --meas-pl ${MEAS_PL} --meas-wl ${MEAS_WL} --meas-pc ${MEAS_PC} \
	      --meas-vt ${MEAS_VT}"

command="./wilson2D ${BASIC_PARAMS} ${HMC_PARAMS} ${EIG_PARAMS} ${MEASUREMENTS}"

echo $command
$command

