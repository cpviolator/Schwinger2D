#!/bin/bash

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum,eig}}

# The RNG seed
SEED=1234

# Verbosities
#-----------------------------------------------
# General
VERBOSITY=1
# CG verbosity
CG_VERBOSITY=0
# Eigensolver verbosity
EIG_VERBOSITY=0

# Physics params
#-----------------------------------------------
# Lattice dims
LX=16
LY=16
# The value of the coupling in the U(1) 2D theory
BETA=4
# Dynamic fermion parameters
# 0 = pure gauge
# 2 = two light degenerate fermions
# 3 = two light degenerate, one heavy fermion
FLAVOURS=0
# Light Fermions (degenerate) mass
MASS=0.05
# Heavy Fermion mass
MASS_HEAVY=0.5

# HMC params
#-----------------------------------------------
# The total number of thermalised HMC iterations to perform.
HMC_ITER=10000
# The number of HMC iterations for thermalisation (accept + accept/reject).
HMC_THERM=100
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=10
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=1000
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=000
# Reverse the gauge fields for ergodicity check
HMC_REVERSE=1000
# HMC time steps in the integration 
HMC_NSTEP=2
# HMC inner time steps in the integration 
HMC_INNER_NSTEP=2
# Degree of polynomial for AlgRemez
HMC_AR_DEGREE=12
# Precision for AlgRemez (GMP)
HMC_AR_GMP_PREC=40
# HMC trajectory time
HMC_TAU=1.5
# Integrator type: leapfrog = 0, FGI = 1
# FYI, for optimal FLOP usage aim for:
# 70% acceptance with Leapfrog
# 90% acceptance with FGI
HMC_INTEGRATOR=1
# Maximum CG iterations
MAX_CG_ITER=1000
# CG residual tolerance 
CG_TOL=1e-10

# Eigensolver parameters
#-----------------------------------------------
# Use the FEAST eigensolver, else IRAM
EIG_USE_FEAST=1
# FEAST M0
EIG_FEAST_M0=40
# FEAST N contours
EIG_FEAST_NCONTOUR=8
# FEAST Max Eval
EIG_FEAST_EMAX=0.8
# FEAST Use previous space as init guess
EIG_FEAST_INIT_GUESS=1

NKR=128
NEV=64
NCONV=32
NDEFL=32

# IRAM MG projector
X_BLK=8
Y_BLK=8
N_LOW=8
# Tolerance on eigenvector residua
EIG_TOL=1e-7
# Maximum restart iterations
MAX_IRAM_ITER=1000
# Operator to solve
EIG_OP=MdagM
# Spectrum to compute
EIG_SPEC=SR
# Compute a deflation space at the start of each
# HMC trajectory and use it throughout the HMC integration
EIG_DEFLATE=1
# Inspect the eigenspectrum at each instance that
# the gauge field is updated in the HMC integrator
EIG_INSPECTION=1
# When deflating, use the compressed space
EIG_USE_COMP_SPACE=false

# Measurements
#-----------------------------------------------
# Smearing algorithm
SMEAR_TYPE=APE
# Number of smearing hits to perform when measuring topology
APE_ITER=5
# The alpha value in the APE smearing
APE_ALPHA=0.6
# Number of smearing hits to perform when measuring topology
WILSON_STEPS=5
# The alpha value in the APE smearing
WILSON_TIME=0.2
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
            --eig-low-modes ${N_LOW} --eig-verbosity ${EIG_VERBOSITY} --eig-deflate ${EIG_DEFLATE} --eig-inspection ${EIG_INSPECTION} \
            --eig-use-comp-space ${EIG_USE_COMP_SPACE} --eig-feast ${EIG_USE_FEAST} \
            --eig-feast-M0 ${EIG_FEAST_M0} --eig-feast-Emax ${EIG_FEAST_EMAX} --eig-feast-Ncontour ${EIG_FEAST_NCONTOUR} --eig-feast-init-guess ${EIG_FEAST_INIT_GUESS} "

MEASUREMENTS="--ape-alpha ${APE_ALPHA} --ape-iter ${APE_ITER} --wilson-steps ${WILSON_STEPS} --wilson-time ${WILSON_TIME} --smear-type ${SMEAR_TYPE} \
	      --meas-pl ${MEAS_PL} --meas-wl ${MEAS_WL} --meas-pc ${MEAS_PC} --meas-vt ${MEAS_VT} "

command="./wilson2D ${BASIC_PARAMS} ${HMC_PARAMS} ${EIG_PARAMS} ${MEASUREMENTS}"

echo $command
$command
