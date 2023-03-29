#!/bin/bash

export OMP_NUM_THREADS=1

mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum,eig}}

# The RNG seed
SEED=1234

# General verbosity
VERBOSITY=0

# Physics params
#-----------------------------------------------
# Lattice dims
LX=8
LY=8
# The value of the coupling in the U(1) 2D theory
BETA=7.0
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
HMC_ITER=100
# The number of HMC iterations for thermalisation (accept + accept/reject).
HMC_THERM=100
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=10
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=100
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=000
# Reverse the gauge fields for ergodicity check
HMC_REVERSE=10000
# HMC time steps in the integration 
HMC_NSTEP=100
# HMC inner time steps in the integration 
HMC_INNER_NSTEP=1
# Degree of polynomial for AlgRemez
HMC_AR_DEGREE=12
# Precision for AlgRemez (GMP)
HMC_AR_GMP_PREC=40
# HMC trajectory time
HMC_TAU=1.0
# Integrator type: leapfrog = 0, omelyan = 1, fgi = 2
# FYI, aim for 70% acceptance with Leapfrog
# and 90% with FGI or OMELYAN for optimal FLOP usage
HMC_INTEGRATOR=0
# Sampler type: HMC = 0, MCHMC = 1
HMC_SAMPLER=1
# beta-eps coefficient
HMC_BETA_EPS=0.2
# Maximum CG iterations
MAX_CG_ITER=10000
# CG residual tolerance 
CG_TOL=1e-9
# CG verbosity
CG_VERBOSITY=false
# Sampler distinction
SAMPLER_DISTINCTION=false

# Eigensolver parameters
#-----------------------------------------------
NKR=128
NEV=128
NCONV=128
NDEFL=128
# IRAM MG projector
X_BLK=4
Y_BLK=4
N_LOW=8
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
# Compute a deflation space at the start of each
# HMC trajectory and use it throughout the HMC integration
EIG_DEFLATE=false
# Inspect the eigenspectrum at each instance that
# the gauge field is updated in the HMC integrator
EIG_INSPECTION=false
# When deflating, use the compressed space
EIG_USE_COMP_SPACE=false

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
            --hmc-integrator ${HMC_INTEGRATOR} --hmc-sampler ${HMC_SAMPLER} --beta-eps ${HMC_BETA_EPS} \
	    --pfe-degree ${HMC_AR_DEGREE} --pfe-prec ${HMC_AR_GMP_PREC} \
	    --cg-max-iter ${MAX_CG_ITER} --cg-tol ${CG_TOL} --cg-verbosity ${CG_VERBOSITY} \
        --sampler-distinction ${SAMPLER_DISTINCTION}"

EIG_PARAMS="--eig-n-ev ${NEV} --eig-n-kr ${NKR} --eig-n-conv ${NCONV} --eig-n-deflate ${NDEFL} --eig-max-restarts ${MAX_IRAM_ITER} \
            --eig-tol ${EIG_TOL} --eig-operator ${EIG_OP} --eig-spectrum ${EIG_SPEC} --eig-block-scheme ${X_BLK} ${Y_BLK} \
            --eig-low-modes ${N_LOW} --eig-verbosity ${EIG_VERBOSITY} --eig-deflate ${EIG_DEFLATE} --eig-inspection ${EIG_INSPECTION} --eig-use-comp-space ${EIG_USE_COMP_SPACE} "

MEASUREMENTS="--ape-alpha ${APE_ALPHA} --ape-iter ${APE_ITER} --meas-pl ${MEAS_PL} --meas-wl ${MEAS_WL} --meas-pc ${MEAS_PC} \
	      --meas-vt ${MEAS_VT}"

command="./wilson2D ${BASIC_PARAMS} ${HMC_PARAMS} ${EIG_PARAMS} ${MEASUREMENTS}"

echo $command
$command

