# Build

Make a build directory and enter it. Then invoke the CMake CL, configure the build, and make. For a quick set up:

```
1. mkdir build
2. cmake /path/to/Schwinger2D/
3. make -j
```

There is a launcher script in the `wilson` directory. Place the launcher script and the executable `wilson2D`
in the same directory in which you wish to run and launch with `./launcher.sh`

# Tutorial

## Getting help

The user may get help on using the executable by running:
```
./wilson2D --help


This is an exhaustive list of run-time options. Those not set at
run time will be given a default value. Important input options will
dumped to stdout at execution. Please ensure that they are sensible!

Default values are given in parentheses

KEY: <N>     = integer
     <float> = decimal
     <bool>  = true/false
GENERAL PARAMS
--help                           Print this message.
--seed <N>                       Sets the RNG seed (1234).
--verbosity <bool>               Sets the verbosity as either verbose or quiet (true).

PHYSICS PARAMS
--beta <float>                   Sets the beta value (3.0).
--dim <N> <N>                    Sets the lattice dimensions (32 32).
--mass <float>                   Sets the light, degenerate quark mass (0.1).
--mass-heavy <float>             Sets the heavy quark mass (0.5).
--flavours <N>                   Sets the number of fermion flavours in the simulation (2).
                                 0=pure gauge, 2=two degenerate, 3=two degenerate + one heavy.

HMC PARAMS
--hmc-traj-length <float>        Sets the HMC trajectory length (1.0).
--hmc-n-step <N>                 Sets the number of HMC outer steps (4).
--hmc-inner-step <N>             Sets the number of HMC inner steps (1).
--hmc-n-trajectories <N>         Sets the number of HMC trajectories after thermalisation (1000).
--hmc-therm <N>                  Sets the number of HMC thermalisation steps (250).
--hmc-checkpoint <N>             Checkpoint the gauge field with this frequency (100).
--hmc-checkpoint-start <N>       Start the HMC from this saved checkpoint (0).
--hmc-skip <N>                   Skip this number of trajectories between measurements (5).
--hmc-reverse <N>                Perform a reversibility check with this frequency (100).
--hmc-integrator <LEAPFROG/FGI>  Sets the HMC integrator (FGI).
--pfe-degree <N>                 Degree of the rational polynomial required for heavy mass fermion (15).
--pfe-prec <N>                   GMP bit-wise precision of the rational polynomal computation (50).
--cg-max-iter <N>                Maximum CG iterations before failure exit (1000).
--cg-tol <float>                 Relative residual norm of CG solution (1e-9).
--cg-verbosity <bool>            Sets CG verbosity as verbose or quiet (false).

EIGENSOLVER PARAMS
--eig-n-ev <N>                   Size of IRAM search space (16).
--eig-n-kr <N>                   Size of IRAM Krylov space (32).
--eig-n-conv <N>                 Number of required converged eigenpairs (16).
--eig-n-deflate <N>              Number of eigenpairs to use in CG deflation. (0)
                                 If set to 0, eigensolver will not be called.
--eig-max-restarts <N>           Maximum number of IRAM restarts (100).
--eig-tol <float>                Residual norm of eigenvectors (1e-9).
--eig-operator <M,Mdag,MdagM, MMdag> Operator to eigensolve (MdagM).
--eig-spectrum <SR, LR, SI, LI, SM, LM> Spectrum to compute (SR)
                                 (S)malest/(L)argest (R)eal/(I)maginary/(M)odulus.
--eig-block-scheme <N> <N>       Size of coarsening scheme for MG projector (2 2).
--eig-low-modes <N>              Number of low eigenmodes used to construct MG projector (16).
--eig-verbosity <bool>           Sets IRAM verbosity as verbose or quiet (false).
--eig-deflate <bool>             Compute a deflation space at the start of the HMC trajectory
                                 and use it throughout the HMC integration (false)
--eig-inspection <bool>          Inspect the eigenspectrum at each call of CG (false).
--eig-use-comp-space <bool>      Use the compressed space for deflation (false).

MEASUREMENT PARAMS
--ape-alpha <float>              Projection coefficient for APE smearing (0.5).
--ape-iter <N>                   Number of APE smearing hits (1).
--meas-pl <bool>                 Measure Polyakov loops every measurement interval (false).
--meas-wl <bool>                 Measure Wilson loops every measurement interval (false).
--meas-pc <bool>                 Measure Pion every measurement interval (false).
--meas-vt <bool>                 Measure Vacuum trace every measurement interval (false).
```

## Simple HMC evolution

The following is example output using the `launcher.sh` script as shipped. In this run,
verbosity is minimal and no eigensolver is called during the HMC evolution.
```
./launcher.sh
./wilson2D --seed 1234 --verbosity 0 --beta 4.0 --dim 8 8 --mass 0.05 --mass-heavy 0.5 --flavours 2 --hmc-traj-length 1.0 --hmc-n-step 4 --hmc-inner-step 1 --hmc-n-trajectories 1000 --hmc-therm 100 --hmc-checkpoint 100 --hmc-checkpoint-start 000 --hmc-skip 10 --hmc-reverse 100 --hmc-integrator 1 --pfe-degree 12 --pfe-prec 40 --cg-max-iter 10000 --cg-tol 1e-9 --cg-verbosity false --eig-n-ev 128 --eig-n-kr 128 --eig-n-conv 128 --eig-n-deflate 128 --eig-max-restarts 1000 --eig-tol 1e-10 --eig-operator MdagM --eig-spectrum SR --eig-block-scheme 4 4 --eig-low-modes 8 --eig-verbosity false --eig-deflate false --eig-inspection false --eig-use-comp-space false --ape-alpha 0.5 --ape-iter 0 --meas-pl 0 --meas-wl 0 --meas-pc 0 --meas-vt 0

********************************
*      Parameter status        *
********************************

Physics:      X dim = 8
              Y dim = 8
              Beta = 4
              Flavours = 2
              Mass = 0.05
HMC:          Therm Sweeps = 100
              Data Points = 1000
              Start Point = 0
              Integrator = FGI
              Trajectory Length = 1
              Trajectory Steps = 4
              Inner Trajectory Steps = 1
Measurements: APE iter = 0
              APE alpha = 0.5
              Polyakov loops = false
              Wilson loops = false
              Pion = false
              Vacuum trace = false
Deflation:    OFF 
1 0.0936080000000000 
2 0.1741220000000000 
3 0.2536700000000000 
4 0.3316350000000000 
5 0.4115390000000000 
6 0.4919980000000000 
7 0.5747850000000000 
8 0.6548130000000000 
9 0.7341130000000000 
10 0.8132340000000000 
11 0.8929580000000000 
12 0.9726020000000000 
13 1.0504549999999999 
14 1.1310070000000001 
15 1.2124100000000000 
16 1.2919160000000001 
17 1.3732000000000000 
18 1.4552419999999999 
19 1.5367800000000000 
20 1.6173649999999999 
21 1.6960409999999999 
22 1.7754090000000000 
23 1.8539930000000000 
24 1.9329769999999999 
25 2.0133749999999999 
26 2.0915919999999999 
27 2.1703869999999998 
28 2.2493850000000002 
29 2.3291789999999999 
30 2.4078349999999999 
31 2.4872470000000000 
32 2.5667140000000002 
33 2.6459100000000002 
34 2.7240199999999999 
35 2.8015189999999999 
36 2.8813390000000001 
37 2.9609990000000002 
38 3.0368800000000000 
39 3.1169769999999999 
40 3.1963080000000001 
41 3.2768259999999998 
42 3.3585950000000002 
43 3.4392839999999998 
44 3.5194179999999999 
45 3.5999140000000001 
46 3.6812870000000002 
47 3.7625709999999999 
48 3.8421780000000001 
49 3.9209339999999999 
50 3.9988280000000000 
51 4.0768639999999996 
52 4.1551950000000000 
53 4.2350339999999997 
54 4.3156020000000002 
55 4.3932599999999997 
56 4.4717690000000001 
57 4.5519090000000002 
58 4.6291019999999996 
59 4.7029870000000003 
60 4.7825270000000000 
61 4.8618860000000002 
62 4.9397760000000002 
63 5.0191270000000001 
64 5.0986330000000004 
65 5.1781050000000004 
66 5.2583940000000000 
67 5.3364050000000001 
68 5.4150660000000004 
69 5.4975690000000004 
70 5.5801959999999999 
71 5.6611289999999999 
72 5.7412869999999998 
73 5.8209000000000000 
74 5.9017200000000001 
75 5.9815090000000000 
76 6.0582370000000001 
77 6.1362310000000004 
78 6.2141320000000002 
79 6.2944060000000004 
80 6.3734050000000000 
81 6.4525210000000000 
82 6.5297780000000003 
83 6.6096870000000001 
84 6.6912409999999998 
85 6.7705979999999997 
86 6.8509989999999998 
87 6.9297440000000003 
88 7.0077749999999996 
89 7.0868830000000003 
90 7.1653339999999996 
91 7.2447800000000004 
92 7.3240650000000000 
93 7.4037090000000001 
94 7.4816010000000004 
95 7.5595929999999996 
96 7.6366280000000000 
97 7.7139239999999996 
98 7.7907570000000002 
99 7.8691060000000004 
100 7.9474030000000004 
IO: Writing gauge gauge/gauge_LX8_LY8_B4.000000_MLight0.050000_traj100.dat
110 8.8594539999999995 0.8623844242142911 0.9750070731089584 1.0557246446900881 1.0264783056822036 -0.0542273981476740 -0.0228728973515985 1.0000000000000000 0
120 9.7687899999999992 0.8791125455313098 0.9187520855636728 0.9948201413012048 1.0249956833030744 +0.0051933206744081 -0.0222970963381044 1.0000000000000000 0
130 10.6789790000000000 0.8789675364546203 0.9015712575669332 1.0049781162461882 0.9958443580266855 -0.0049657663945482 -0.0240599424880382 0.9666666666666667 0
140 11.6012529999999998 0.8162623211218231 0.8883692947261455 0.8953758535109170 0.9928087217589429 +0.1105117008038974 -0.0132901653819367 0.9750000000000000 1
150 12.5483010000000004 0.8273569505172855 0.8798395828094243 0.9878624797134208 1.0591579805584892 +0.0122117814956368 -0.0384775611208744 0.9800000000000000 1
160 13.5007310000000000 0.8623661958862654 0.8772137130401718 0.6634393890614280 0.9591646157559440 +0.4103177797044566 0.0022415282700313 0.9166666666666666 1
170 14.4423220000000008 0.8755492879225648 0.8757882363324826 1.1770177752828477 0.9603581200520622 -0.1629839303585356 0.0095986509399365 0.9285714285714286 1
180 15.3837150000000005 0.8488887400889705 0.8751895516145295 0.5047065630010715 0.9681904753195860 +0.6837780819537329 0.0054467788774112 0.9250000000000000 1
190 16.3373139999999992 0.8518511853003609 0.8732416724640117 0.9596082375647987 0.9318024788683407 +0.0412301636799555 0.0121711238576836 0.9000000000000000 1
H0 = 229.8306817994547089
H1 = 229.8403330358098628
H2 = 229.8306818017462376
H2 - H0 = 2.2915287445357535e-09
(H2 - H0)/H0 percent error = 9.970508404684e-10
L2 norm of gauge(t0) - evolve_gauge(t0->t1->t0) = 1.9062085890424138e-21
L2 norm of gauge(t0) - flowed_gauge(t0->t1->t0) = 9.1893321174617302e-09
IO: Writing gauge gauge/gauge_LX8_LY8_B4.000000_MLight0.050000_traj200.dat
200 17.4598230000000001 0.8836837405677062 0.8717285863450196 0.9922978283872238 0.9439015124795072 +0.0077319865282846 0.0063276318842207 0.9100000000000000 0
...

```
Points to note:
1. The command is printed to the terminal
2. The physics and algorithmic parameters are printed.
3. 100 thermalisation sweeps are performed and always accepted.
4. After the thermalisation sweeps are performed, `--hmc-skip` sweeps are
   performed and then measurements are taken.
5. After thermalisation, from left to right, the data in the terminal output is
   (HMC iteration) (wall clock) (lattice action) (average lattice action) (exp(-dH)) (dH) (average dH) (Average acceptance) (topological charge)
6. The gauge field is saved every `--hmc-checkpoint` iteratiions of HMC.
7. A reversibility check is performed every `--hmc-reverse` iterations of HMC to ensure ergodicity.

## HMC evolution with a deflation space

The following is example output using the `launcher.sh` script with deflation turned on. To do this, we set `--eig-deflate true`.
In this run, verbosity is on and an eigensolver is called at the start the HMC integration to construct a deflation space.
by setting `--hmc-therm 0` and `--hmc-checkpoint-start 100` we pick up the already thermalised lattice
from HMC iteration 100.

```
./launcher.sh
./wilson2D --seed 1234 --verbosity 1 --beta 4.0 --dim 8 8 --mass 0.05 --mass-heavy 0.5 --flavours 2 --hmc-traj-length 1.0 --hmc-n-step 4 --hmc-inner-step 1 --hmc-n-trajectories 1000 --hmc-therm 000 --hmc-checkpoint 100 --hmc-checkpoint-start 100 --hmc-skip 10 --hmc-reverse 100 --hmc-integrator 1 --pfe-degree 12 --pfe-prec 40 --cg-max-iter 10000 --cg-tol 1e-9 --cg-verbosity false --eig-n-ev 128 --eig-n-kr 128 --eig-n-conv 128 --eig-n-deflate 128 --eig-max-restarts 1000 --eig-tol 1e-10 --eig-operator MdagM --eig-spectrum SR --eig-block-scheme 4 4 --eig-low-modes 8 --eig-verbosity false --eig-deflate true --eig-inspection false --eig-use-comp-space false --ape-alpha 0.5 --ape-iter 0 --meas-pl 0 --meas-wl 0 --meas-pc 0 --meas-vt 0

********************************
*      Parameter status        *
********************************

Physics:      X dim = 8
              Y dim = 8
              Beta = 4
              Flavours = 2
              Mass = 0.05
HMC:          Therm Sweeps = 0
              Data Points = 1000
              Start Point = 100
              Integrator = FGI
              Trajectory Length = 1
              Trajectory Steps = 4
              Inner Trajectory Steps = 1
Measurements: APE iter = 0
              APE alpha = 0.5
              Polyakov loops = false
              Wilson loops = false
              Pion = false
              Vacuum trace = false
Deflation:    nkv = 128
              nev = 128
              ndefl = 128
              tol = 1e-10
              maxiter = 1000
IO: Reading gauge gauge/gauge_LX8_LY8_B4.000000_MLight0.050000_traj100.dat
CG: Switching off eigensolver
CG: Converged iter = 73, res = 1.8854125564414205e-08, truersq = 7.1365884137920052e-10
CG: Restoring eigensolver
CG: Computing deflation space for CG
IRAM: 0128 converged eigenvalues at iter 0
IRAM: Computed the requested 128 vectors with a 128 search space and a 128 Krylov space in 1 restart_steps and 128 OPs in 3.034590e+00 secs.
IRAM: EigValue[0000]: ||(+1.47247383e-01, -1.27097294e-18)|| = +1.47247383e-01 residual 2.29569921e-14
IRAM: EigValue[0001]: ||(+1.50830439e-01, -7.64913103e-18)|| = +1.50830439e-01 residual 2.69911376e-14
IRAM: EigValue[0002]: ||(+2.11832954e-01, -1.41793315e-18)|| = +2.11832954e-01 residual 2.57323302e-14
IRAM: EigValue[0003]: ||(+2.45975362e-01, -4.05305265e-19)|| = +2.45975362e-01 residual 2.53617511e-14
IRAM: EigValue[0004]: ||(+3.68019589e-01, -3.40507245e-19)|| = +3.68019589e-01 residual 2.14851256e-14
IRAM: EigValue[0005]: ||(+3.92991823e-01, -6.12658931e-18)|| = +3.92991823e-01 residual 2.33811217e-14
IRAM: EigValue[0006]: ||(+5.54948852e-01, -3.91498628e-18)|| = +5.54948852e-01 residual 2.13155557e-14
IRAM: EigValue[0007]: ||(+5.65607264e-01, -1.13317987e-17)|| = +5.65607264e-01 residual 1.98934620e-14
IRAM: EigValue[0008]: ||(+1.04912932e+00, +2.35271871e-17)|| = +1.04912932e+00 residual 1.68023792e-14
IRAM: EigValue[0009]: ||(+1.07001379e+00, -1.67153482e-17)|| = +1.07001379e+00 residual 2.01038351e-14
IRAM: EigValue[0010]: ||(+1.11552921e+00, -1.23204648e-17)|| = +1.11552921e+00 residual 1.90504169e-14
IRAM: EigValue[0011]: ||(+1.14993109e+00, +3.49146981e-18)|| = +1.14993109e+00 residual 1.96571553e-14
IRAM: EigValue[0012]: ||(+1.23041416e+00, -1.69237183e-17)|| = +1.23041416e+00 residual 1.80378770e-14
IRAM: EigValue[0013]: ||(+1.29259816e+00, +1.50250939e-17)|| = +1.29259816e+00 residual 1.54578848e-14
IRAM: EigValue[0014]: ||(+1.41098354e+00, +1.26783892e-17)|| = +1.41098354e+00 residual 1.98222484e-14
IRAM: EigValue[0015]: ||(+1.42343230e+00, -1.55108673e-17)|| = +1.42343230e+00 residual 1.70445069e-14
IRAM: EigValue[0016]: ||(+1.56346758e+00, -1.55938766e-18)|| = +1.56346758e+00 residual 1.48283793e-14
IRAM: EigValue[0017]: ||(+1.60882700e+00, -5.23127548e-17)|| = +1.60882700e+00 residual 1.31660862e-14
IRAM: EigValue[0018]: ||(+1.76395268e+00, -1.48614048e-17)|| = +1.76395268e+00 residual 1.58078839e-14
IRAM: EigValue[0019]: ||(+1.76531521e+00, -3.22471796e-17)|| = +1.76531521e+00 residual 1.39746131e-14
IRAM: EigValue[0020]: ||(+1.87848659e+00, -2.02576400e-17)|| = +1.87848659e+00 residual 1.67691019e-14
IRAM: EigValue[0021]: ||(+1.92404791e+00, +7.41390998e-17)|| = +1.92404791e+00 residual 1.38480001e-14
IRAM: EigValue[0022]: ||(+2.06186721e+00, -2.85483985e-17)|| = +2.06186721e+00 residual 1.44940463e-14
IRAM: EigValue[0023]: ||(+2.11825709e+00, -1.81052234e-17)|| = +2.11825709e+00 residual 1.24261723e-14
IRAM: EigValue[0024]: ||(+2.43526104e+00, +5.70730800e-18)|| = +2.43526104e+00 residual 1.22378969e-14
IRAM: EigValue[0025]: ||(+2.48301514e+00, -3.62486162e-17)|| = +2.48301514e+00 residual 1.12697546e-14
IRAM: EigValue[0026]: ||(+2.65203822e+00, -8.39409651e-17)|| = +2.65203822e+00 residual 8.85491000e-15
IRAM: EigValue[0027]: ||(+2.66989848e+00, +5.09306353e-18)|| = +2.66989848e+00 residual 9.49306114e-15
IRAM: EigValue[0028]: ||(+2.77779529e+00, -2.33594746e-17)|| = +2.77779529e+00 residual 8.50067834e-15
IRAM: EigValue[0029]: ||(+2.77955875e+00, -4.86148207e-17)|| = +2.77955875e+00 residual 7.95849447e-15
IRAM: EigValue[0030]: ||(+2.82578441e+00, -1.71922277e-17)|| = +2.82578441e+00 residual 1.14363921e-14
IRAM: EigValue[0031]: ||(+2.86653354e+00, +3.33087236e-17)|| = +2.86653354e+00 residual 9.50896173e-15
IRAM: EigValue[0032]: ||(+3.04239484e+00, +1.43267153e-17)|| = +3.04239484e+00 residual 9.79983754e-15
IRAM: EigValue[0033]: ||(+3.05268486e+00, -2.06413459e-17)|| = +3.05268486e+00 residual 1.06854806e-14
IRAM: EigValue[0034]: ||(+3.17106693e+00, +4.43209990e-18)|| = +3.17106693e+00 residual 8.83373993e-15
IRAM: EigValue[0035]: ||(+3.18113584e+00, -1.48366291e-17)|| = +3.18113584e+00 residual 7.20299736e-15
IRAM: EigValue[0036]: ||(+3.28758020e+00, +7.69156738e-17)|| = +3.28758020e+00 residual 7.09464092e-15
IRAM: EigValue[0037]: ||(+3.31300072e+00, -4.44624535e-17)|| = +3.31300072e+00 residual 5.35029014e-15
IRAM: EigValue[0038]: ||(+3.39304765e+00, -1.59428541e-17)|| = +3.39304765e+00 residual 8.27434988e-15
IRAM: EigValue[0039]: ||(+3.42404793e+00, +2.67548459e-17)|| = +3.42404793e+00 residual 6.43378580e-15
IRAM: EigValue[0040]: ||(+3.51260753e+00, -1.31525158e-17)|| = +3.51260753e+00 residual 6.86596298e-15
IRAM: EigValue[0041]: ||(+3.53421851e+00, -3.98444298e-18)|| = +3.53421851e+00 residual 6.15594653e-15
IRAM: EigValue[0042]: ||(+3.62721271e+00, +1.71237239e-17)|| = +3.62721271e+00 residual 5.54739556e-15
IRAM: EigValue[0043]: ||(+3.66034520e+00, +4.16774091e-17)|| = +3.66034520e+00 residual 5.41212897e-15
IRAM: EigValue[0044]: ||(+3.78870735e+00, -2.51196091e-17)|| = +3.78870735e+00 residual 3.71915283e-15
IRAM: EigValue[0045]: ||(+3.79658510e+00, -5.67130910e-18)|| = +3.79658510e+00 residual 7.33039657e-15
IRAM: EigValue[0046]: ||(+3.91003502e+00, +3.66028347e-17)|| = +3.91003502e+00 residual 1.07821184e-14
IRAM: EigValue[0047]: ||(+3.98064995e+00, +4.01866311e-17)|| = +3.98064995e+00 residual 1.14912841e-14
IRAM: EigValue[0048]: ||(+4.09648893e+00, -7.20825038e-19)|| = +4.09648893e+00 residual 8.21555224e-15
IRAM: EigValue[0049]: ||(+4.12511494e+00, -6.59669259e-18)|| = +4.12511494e+00 residual 3.72243994e-15
IRAM: EigValue[0050]: ||(+4.21934741e+00, -2.18305801e-17)|| = +4.21934741e+00 residual 3.49597382e-15
IRAM: EigValue[0051]: ||(+4.23591669e+00, +5.10921803e-17)|| = +4.23591669e+00 residual 3.98287535e-15
IRAM: EigValue[0052]: ||(+4.35205969e+00, -1.35626916e-17)|| = +4.35205969e+00 residual 5.95725624e-15
IRAM: EigValue[0053]: ||(+4.36940878e+00, +2.93776437e-17)|| = +4.36940878e+00 residual 5.93283472e-15
IRAM: EigValue[0054]: ||(+4.47983461e+00, +6.25813352e-17)|| = +4.47983461e+00 residual 9.05394807e-15
IRAM: EigValue[0055]: ||(+4.50794453e+00, +2.42209071e-17)|| = +4.50794453e+00 residual 5.85697288e-15
IRAM: EigValue[0056]: ||(+4.64680045e+00, +8.90698554e-17)|| = +4.64680045e+00 residual 7.21365150e-15
IRAM: EigValue[0057]: ||(+4.66233520e+00, -1.21600050e-17)|| = +4.66233520e+00 residual 8.58332807e-15
IRAM: EigValue[0058]: ||(+4.74084958e+00, -4.66037528e-18)|| = +4.74084958e+00 residual 6.07351387e-15
IRAM: EigValue[0059]: ||(+4.77252548e+00, +5.65377552e-17)|| = +4.77252548e+00 residual 6.03233303e-15
IRAM: EigValue[0060]: ||(+4.89822464e+00, -1.43079111e-16)|| = +4.89822464e+00 residual 6.41584905e-15
IRAM: EigValue[0061]: ||(+4.90575414e+00, -3.04940331e-17)|| = +4.90575414e+00 residual 5.82598134e-15
IRAM: EigValue[0062]: ||(+5.11946117e+00, +3.63715948e-17)|| = +5.11946117e+00 residual 8.27546771e-15
IRAM: EigValue[0063]: ||(+5.14192476e+00, -2.58640981e-17)|| = +5.14192476e+00 residual 1.66863039e-14
IRAM: EigValue[0064]: ||(+5.26808093e+00, -4.77455532e-17)|| = +5.26808093e+00 residual 1.42973794e-14
IRAM: EigValue[0065]: ||(+5.27037669e+00, +5.11078505e-17)|| = +5.27037669e+00 residual 2.31659814e-14
IRAM: EigValue[0066]: ||(+5.40338032e+00, -4.29447822e-17)|| = +5.40338032e+00 residual 6.75500097e-15
IRAM: EigValue[0067]: ||(+5.40880935e+00, -4.71983699e-17)|| = +5.40880935e+00 residual 7.53971003e-15
IRAM: EigValue[0068]: ||(+5.48192030e+00, +1.35609975e-17)|| = +5.48192030e+00 residual 6.23312466e-15
IRAM: EigValue[0069]: ||(+5.50733702e+00, +5.79724172e-17)|| = +5.50733702e+00 residual 7.41829990e-15
IRAM: EigValue[0070]: ||(+5.79416201e+00, +2.31578808e-17)|| = +5.79416201e+00 residual 8.96932527e-15
IRAM: EigValue[0071]: ||(+5.80781865e+00, +1.96901279e-17)|| = +5.80781865e+00 residual 1.09210505e-14
IRAM: EigValue[0072]: ||(+6.08130698e+00, +5.79573824e-17)|| = +6.08130698e+00 residual 1.44219784e-14
IRAM: EigValue[0073]: ||(+6.15590477e+00, -4.92281255e-17)|| = +6.15590477e+00 residual 1.43469160e-14
IRAM: EigValue[0074]: ||(+6.29049999e+00, -4.35984799e-17)|| = +6.29049999e+00 residual 1.20961277e-14
IRAM: EigValue[0075]: ||(+6.33168378e+00, -1.91632734e-17)|| = +6.33168378e+00 residual 1.10554097e-14
IRAM: EigValue[0076]: ||(+6.46134095e+00, +6.11434968e-17)|| = +6.46134095e+00 residual 1.99101570e-14
IRAM: EigValue[0077]: ||(+6.47683109e+00, +3.34747421e-17)|| = +6.47683109e+00 residual 1.05074508e-14
IRAM: EigValue[0078]: ||(+6.70307632e+00, +1.53507781e-17)|| = +6.70307632e+00 residual 1.14755168e-14
IRAM: EigValue[0079]: ||(+6.73923841e+00, +7.88037102e-18)|| = +6.73923841e+00 residual 2.03711332e-14
IRAM: EigValue[0080]: ||(+6.89014735e+00, -1.75793218e-17)|| = +6.89014735e+00 residual 1.56354396e-14
IRAM: EigValue[0081]: ||(+6.96401323e+00, -5.53832493e-18)|| = +6.96401323e+00 residual 1.74853751e-14
IRAM: EigValue[0082]: ||(+7.20994999e+00, -7.58492593e-17)|| = +7.20994999e+00 residual 2.91127779e-14
IRAM: EigValue[0083]: ||(+7.21831072e+00, +1.24626075e-16)|| = +7.21831072e+00 residual 2.15661683e-14
IRAM: EigValue[0084]: ||(+7.41295610e+00, +6.64751457e-18)|| = +7.41295610e+00 residual 1.74381735e-14
IRAM: EigValue[0085]: ||(+7.47203720e+00, -7.04375658e-17)|| = +7.47203720e+00 residual 2.65591265e-14
IRAM: EigValue[0086]: ||(+7.68588751e+00, +1.38204860e-16)|| = +7.68588751e+00 residual 1.62298602e-14
IRAM: EigValue[0087]: ||(+7.69280102e+00, -1.99832013e-17)|| = +7.69280102e+00 residual 1.44666247e-14
IRAM: EigValue[0088]: ||(+7.97301663e+00, -4.40999234e-17)|| = +7.97301663e+00 residual 2.89745606e-14
IRAM: EigValue[0089]: ||(+8.02465232e+00, -6.17448902e-17)|| = +8.02465232e+00 residual 2.75854911e-14
IRAM: EigValue[0090]: ||(+8.24103952e+00, -2.01198383e-16)|| = +8.24103952e+00 residual 3.14525472e-14
IRAM: EigValue[0091]: ||(+8.26625360e+00, +4.15012263e-17)|| = +8.26625360e+00 residual 3.71165109e-14
IRAM: EigValue[0092]: ||(+8.58205939e+00, +1.31147805e-16)|| = +8.58205939e+00 residual 2.10978271e-14
IRAM: EigValue[0093]: ||(+8.58383048e+00, +4.38351197e-17)|| = +8.58383048e+00 residual 2.10583335e-14
IRAM: EigValue[0094]: ||(+8.77455568e+00, +1.46311813e-16)|| = +8.77455568e+00 residual 4.50516386e-14
IRAM: EigValue[0095]: ||(+8.78283446e+00, -1.31056114e-17)|| = +8.78283446e+00 residual 2.65419384e-14
IRAM: EigValue[0096]: ||(+8.94542528e+00, +2.43945489e-19)|| = +8.94542528e+00 residual 3.47067841e-14
IRAM: EigValue[0097]: ||(+8.96797202e+00, -7.49624158e-18)|| = +8.96797202e+00 residual 2.51388442e-14
IRAM: EigValue[0098]: ||(+9.46488406e+00, +1.77843038e-17)|| = +9.46488406e+00 residual 2.76027835e-14
IRAM: EigValue[0099]: ||(+9.58532672e+00, -1.47715770e-16)|| = +9.58532672e+00 residual 1.73207304e-14
IRAM: EigValue[0100]: ||(+9.85161258e+00, -8.58465774e-17)|| = +9.85161258e+00 residual 3.99177786e-14
IRAM: EigValue[0101]: ||(+9.94591340e+00, -1.15311677e-16)|| = +9.94591340e+00 residual 2.52192855e-14
IRAM: EigValue[0102]: ||(+1.04019207e+01, -3.50796578e-17)|| = +1.04019207e+01 residual 5.41354481e-14
IRAM: EigValue[0103]: ||(+1.04305851e+01, +2.85819621e-17)|| = +1.04305851e+01 residual 5.96901454e-14
IRAM: EigValue[0104]: ||(+1.10271537e+01, -2.24395968e-16)|| = +1.10271537e+01 residual 5.37708105e-14
IRAM: EigValue[0105]: ||(+1.10950302e+01, -7.83111606e-17)|| = +1.10950302e+01 residual 6.35357350e-14
IRAM: EigValue[0106]: ||(+1.13732913e+01, -1.60498344e-16)|| = +1.13732913e+01 residual 7.64579835e-14
IRAM: EigValue[0107]: ||(+1.14375919e+01, -4.25325948e-17)|| = +1.14375919e+01 residual 4.59417278e-14
IRAM: EigValue[0108]: ||(+1.16845873e+01, +5.03798256e-17)|| = +1.16845873e+01 residual 3.21106950e-14
IRAM: EigValue[0109]: ||(+1.16894539e+01, +4.08005712e-17)|| = +1.16894539e+01 residual 3.84771903e-14
IRAM: EigValue[0110]: ||(+1.21169663e+01, -1.15149047e-16)|| = +1.21169663e+01 residual 5.24077081e-14
IRAM: EigValue[0111]: ||(+1.21698762e+01, -9.79949357e-17)|| = +1.21698762e+01 residual 4.52715090e-14
IRAM: EigValue[0112]: ||(+1.24142685e+01, -6.69274613e-17)|| = +1.24142685e+01 residual 6.03358411e-14
IRAM: EigValue[0113]: ||(+1.24906008e+01, +1.06453407e-16)|| = +1.24906008e+01 residual 3.52001770e-14
IRAM: EigValue[0114]: ||(+1.28211693e+01, +3.31520225e-17)|| = +1.28211693e+01 residual 5.70569970e-14
IRAM: EigValue[0115]: ||(+1.28789109e+01, -9.56605129e-17)|| = +1.28789109e+01 residual 8.13617633e-14
IRAM: EigValue[0116]: ||(+1.30176504e+01, -2.38304249e-17)|| = +1.30176504e+01 residual 8.04161166e-14
IRAM: EigValue[0117]: ||(+1.30594773e+01, +2.07950400e-16)|| = +1.30594773e+01 residual 6.24653707e-14
IRAM: EigValue[0118]: ||(+1.32423153e+01, +1.32865588e-16)|| = +1.32423153e+01 residual 3.72623627e-14
IRAM: EigValue[0119]: ||(+1.32507605e+01, -1.50839627e-17)|| = +1.32507605e+01 residual 2.90334043e-14
IRAM: EigValue[0120]: ||(+1.47553315e+01, -1.37836877e-17)|| = +1.47553315e+01 residual 7.85481785e-14
IRAM: EigValue[0121]: ||(+1.47720854e+01, -2.59395370e-17)|| = +1.47720854e+01 residual 4.17589089e-14
IRAM: EigValue[0122]: ||(+1.50895773e+01, +4.39876915e-17)|| = +1.50895773e+01 residual 9.07199857e-14
IRAM: EigValue[0123]: ||(+1.51572583e+01, -3.24447500e-17)|| = +1.51572583e+01 residual 5.97064640e-14
IRAM: EigValue[0124]: ||(+1.54360435e+01, +1.59515786e-16)|| = +1.54360435e+01 residual 7.54230866e-14
IRAM: EigValue[0125]: ||(+1.54604865e+01, +1.29416470e-16)|| = +1.54604865e+01 residual 8.11739476e-14
IRAM: EigValue[0126]: ||(+1.57060762e+01, -2.29689924e-17)|| = +1.57060762e+01 residual 1.41347597e-13
IRAM: EigValue[0127]: ||(+1.57213820e+01, -1.14671320e-16)|| = +1.57213820e+01 residual 5.23964257e-14
IRAM: Post Compression EigValue[0000]: ||(+1.47247383e-01, -5.62006361e-19)|| = +1.47247383e-01 residual 2.30063126e-14
IRAM: Post Compression EigValue[0001]: ||(+1.50830439e-01, -1.13576530e-17)|| = +1.50830439e-01 residual 2.68872021e-14
IRAM: Post Compression EigValue[0002]: ||(+2.11832954e-01, -1.15238832e-18)|| = +2.11832954e-01 residual 2.57938327e-14
IRAM: Post Compression EigValue[0003]: ||(+2.45975362e-01, -7.37977455e-19)|| = +2.45975362e-01 residual 2.54928873e-14
IRAM: Post Compression EigValue[0004]: ||(+3.68019589e-01, -7.72663454e-18)|| = +3.68019589e-01 residual 2.20021430e-14
IRAM: Post Compression EigValue[0005]: ||(+3.92991823e-01, -4.69933879e-18)|| = +3.92991823e-01 residual 2.38333009e-14
IRAM: Post Compression EigValue[0006]: ||(+5.54948852e-01, -6.33919458e-18)|| = +5.54948852e-01 residual 2.31906234e-14
IRAM: Post Compression EigValue[0007]: ||(+5.65607264e-01, +3.89578246e-18)|| = +5.65607264e-01 residual 2.78135484e-14
IRAM: Post Compression EigValue[0008]: ||(+1.04912932e+00, -1.29596041e-19)|| = +1.04912932e+00 residual 1.97060568e-14
IRAM: Post Compression EigValue[0009]: ||(+1.07001379e+00, +3.38609891e-17)|| = +1.07001379e+00 residual 2.35938360e-14
IRAM: Post Compression EigValue[0010]: ||(+1.11552921e+00, -1.19247945e-17)|| = +1.11552921e+00 residual 2.21176101e-14
IRAM: Post Compression EigValue[0011]: ||(+1.14993109e+00, -1.67475354e-17)|| = +1.14993109e+00 residual 2.31426251e-14
IRAM: Post Compression EigValue[0012]: ||(+1.23041416e+00, -1.60597447e-17)|| = +1.23041416e+00 residual 2.14795964e-14
IRAM: Post Compression EigValue[0013]: ||(+1.29259816e+00, +9.95221361e-18)|| = +1.29259816e+00 residual 2.01565442e-14
IRAM: Post Compression EigValue[0014]: ||(+1.41098354e+00, +3.40303957e-17)|| = +1.41098354e+00 residual 2.31212971e-14
IRAM: Post Compression EigValue[0015]: ||(+1.42343230e+00, -2.10402984e-17)|| = +1.42343230e+00 residual 2.05343705e-14
IRAM: Post Compression EigValue[0016]: ||(+1.56346758e+00, -1.03922472e-17)|| = +1.56346758e+00 residual 1.83409245e-14
IRAM: Post Compression EigValue[0017]: ||(+1.60882700e+00, -2.75624521e-17)|| = +1.60882700e+00 residual 2.64238243e-14
IRAM: Post Compression EigValue[0018]: ||(+1.76395268e+00, +2.62631036e-17)|| = +1.76395268e+00 residual 1.84906394e-14
IRAM: Post Compression EigValue[0019]: ||(+1.76531521e+00, -1.90832288e-17)|| = +1.76531521e+00 residual 2.61345064e-14
IRAM: Post Compression EigValue[0020]: ||(+1.87848659e+00, +1.55921825e-17)|| = +1.87848659e+00 residual 2.06794440e-14
IRAM: Post Compression EigValue[0021]: ||(+1.92404791e+00, -2.25683458e-17)|| = +1.92404791e+00 residual 3.04756392e-14
IRAM: Post Compression EigValue[0022]: ||(+2.06186721e+00, +1.94055248e-18)|| = +2.06186721e+00 residual 2.74781483e-14
IRAM: Post Compression EigValue[0023]: ||(+2.11825709e+00, +1.39225747e-17)|| = +2.11825709e+00 residual 1.59743980e-14
IRAM: Post Compression EigValue[0024]: ||(+2.43526104e+00, -4.60328526e-17)|| = +2.43526104e+00 residual 1.67736859e-14
IRAM: Post Compression EigValue[0025]: ||(+2.48301514e+00, +1.94206655e-17)|| = +2.48301514e+00 residual 2.86943030e-14
IRAM: Post Compression EigValue[0026]: ||(+2.65203822e+00, +3.57735895e-17)|| = +2.65203822e+00 residual 1.82361811e-14
IRAM: Post Compression EigValue[0027]: ||(+2.66989848e+00, -1.63776004e-17)|| = +2.66989848e+00 residual 1.41827750e-14
IRAM: Post Compression EigValue[0028]: ||(+2.77779529e+00, +9.91994166e-17)|| = +2.77779529e+00 residual 1.94811412e-14
IRAM: Post Compression EigValue[0029]: ||(+2.77955875e+00, +1.01806478e-17)|| = +2.77955875e+00 residual 3.05776378e-14
IRAM: Post Compression EigValue[0030]: ||(+2.82578441e+00, -3.76988954e-17)|| = +2.82578441e+00 residual 1.99910724e-14
IRAM: Post Compression EigValue[0031]: ||(+2.86653354e+00, -1.08111897e-16)|| = +2.86653354e+00 residual 2.60205824e-14
IRAM: Post Compression EigValue[0032]: ||(+3.04239484e+00, +2.63325603e-17)|| = +3.04239484e+00 residual 2.92488811e-14
IRAM: Post Compression EigValue[0033]: ||(+3.05268486e+00, +3.88678009e-17)|| = +3.05268486e+00 residual 2.04529134e-14
IRAM: Post Compression EigValue[0034]: ||(+3.17106693e+00, -2.47786783e-17)|| = +3.17106693e+00 residual 1.62466809e-14
IRAM: Post Compression EigValue[0035]: ||(+3.18113584e+00, -5.38628251e-18)|| = +3.18113584e+00 residual 2.85021004e-14
IRAM: Post Compression EigValue[0036]: ||(+3.28758020e+00, +3.70119517e-17)|| = +3.28758020e+00 residual 1.87691382e-14
IRAM: Post Compression EigValue[0037]: ||(+3.31300072e+00, -1.94851459e-17)|| = +3.31300072e+00 residual 2.81181025e-14
IRAM: Post Compression EigValue[0038]: ||(+3.39304765e+00, +6.31039546e-17)|| = +3.39304765e+00 residual 1.89003193e-14
IRAM: Post Compression EigValue[0039]: ||(+3.42404793e+00, +1.45797611e-16)|| = +3.42404793e+00 residual 2.60770739e-14
IRAM: Post Compression EigValue[0040]: ||(+3.51260753e+00, -3.17531476e-17)|| = +3.51260753e+00 residual 1.83548379e-14
IRAM: Post Compression EigValue[0041]: ||(+3.53421851e+00, -5.78235512e-17)|| = +3.53421851e+00 residual 2.53591335e-14
IRAM: Post Compression EigValue[0042]: ||(+3.62721271e+00, -2.70523265e-17)|| = +3.62721271e+00 residual 1.56233200e-14
IRAM: Post Compression EigValue[0043]: ||(+3.66034520e+00, +5.25262071e-17)|| = +3.66034520e+00 residual 2.63435197e-14
IRAM: Post Compression EigValue[0044]: ||(+3.78870735e+00, +6.40492433e-17)|| = +3.78870735e+00 residual 1.65905641e-14
IRAM: Post Compression EigValue[0045]: ||(+3.79658510e+00, -8.15893899e-17)|| = +3.79658510e+00 residual 3.20357293e-14
IRAM: Post Compression EigValue[0046]: ||(+3.91003502e+00, +9.19403442e-17)|| = +3.91003502e+00 residual 2.93010187e-14
IRAM: Post Compression EigValue[0047]: ||(+3.98064995e+00, -2.07912707e-17)|| = +3.98064995e+00 residual 2.09210013e-14
IRAM: Post Compression EigValue[0048]: ||(+4.09648893e+00, +2.35009291e-18)|| = +4.09648893e+00 residual 2.79612416e-14
IRAM: Post Compression EigValue[0049]: ||(+4.12511494e+00, +7.50183200e-17)|| = +4.12511494e+00 residual 1.63820492e-14
IRAM: Post Compression EigValue[0050]: ||(+4.21934741e+00, -8.02724654e-17)|| = +4.21934741e+00 residual 3.12495109e-14
IRAM: Post Compression EigValue[0051]: ||(+4.23591669e+00, -6.12574227e-18)|| = +4.23591669e+00 residual 1.74956624e-14
IRAM: Post Compression EigValue[0052]: ||(+4.35205969e+00, +3.03280147e-17)|| = +4.35205969e+00 residual 2.09325942e-14
IRAM: Post Compression EigValue[0053]: ||(+4.36940878e+00, +1.74725956e-17)|| = +4.36940878e+00 residual 3.66259805e-14
IRAM: Post Compression EigValue[0054]: ||(+4.47983461e+00, +5.21026907e-17)|| = +4.47983461e+00 residual 1.87732452e-14
IRAM: Post Compression EigValue[0055]: ||(+4.50794453e+00, -2.22684962e-17)|| = +4.50794453e+00 residual 2.35994469e-14
IRAM: Post Compression EigValue[0056]: ||(+4.64680045e+00, -9.94557499e-17)|| = +4.64680045e+00 residual 2.02965162e-14
IRAM: Post Compression EigValue[0057]: ||(+4.66233520e+00, -2.58480574e-17)|| = +4.66233520e+00 residual 2.96853997e-14
IRAM: Post Compression EigValue[0058]: ||(+4.74084958e+00, -5.04535175e-17)|| = +4.74084958e+00 residual 3.22077270e-14
IRAM: Post Compression EigValue[0059]: ||(+4.77252548e+00, -6.73458956e-17)|| = +4.77252548e+00 residual 1.98111349e-14
IRAM: Post Compression EigValue[0060]: ||(+4.89822464e+00, +1.61520713e-17)|| = +4.89822464e+00 residual 1.72604452e-14
IRAM: Post Compression EigValue[0061]: ||(+4.90575414e+00, -5.24830084e-17)|| = +4.90575414e+00 residual 4.07955950e-14
IRAM: Post Compression EigValue[0062]: ||(+5.11946117e+00, +2.72304152e-17)|| = +5.11946117e+00 residual 2.30899378e-14
IRAM: Post Compression EigValue[0063]: ||(+5.14192476e+00, +2.23855985e-17)|| = +5.14192476e+00 residual 3.20237862e-14
IRAM: Post Compression EigValue[0064]: ||(+5.26808093e+00, +3.73202717e-18)|| = +5.26808093e+00 residual 2.75186241e-14
IRAM: Post Compression EigValue[0065]: ||(+5.27037669e+00, -1.67151364e-17)|| = +5.27037669e+00 residual 2.44532116e-14
IRAM: Post Compression EigValue[0066]: ||(+5.40338032e+00, -1.85924791e-17)|| = +5.40338032e+00 residual 2.02113252e-14
IRAM: Post Compression EigValue[0067]: ||(+5.40880935e+00, -1.68823831e-16)|| = +5.40880935e+00 residual 3.59915866e-14
IRAM: Post Compression EigValue[0068]: ||(+5.48192030e+00, -1.79198290e-17)|| = +5.48192030e+00 residual 3.02333884e-14
IRAM: Post Compression EigValue[0069]: ||(+5.50733702e+00, -7.55885849e-17)|| = +5.50733702e+00 residual 1.77308892e-14
IRAM: Post Compression EigValue[0070]: ||(+5.79416201e+00, -1.06133228e-17)|| = +5.79416201e+00 residual 1.93065464e-14
IRAM: Post Compression EigValue[0071]: ||(+5.80781865e+00, +1.09013140e-16)|| = +5.80781865e+00 residual 2.68797977e-14
IRAM: Post Compression EigValue[0072]: ||(+6.08130698e+00, -1.21268013e-16)|| = +6.08130698e+00 residual 2.75579805e-14
IRAM: Post Compression EigValue[0073]: ||(+6.15590477e+00, -5.38466256e-17)|| = +6.15590477e+00 residual 3.25354326e-14
IRAM: Post Compression EigValue[0074]: ||(+6.29049999e+00, +5.61989420e-17)|| = +6.29049999e+00 residual 2.55435973e-14
IRAM: Post Compression EigValue[0075]: ||(+6.33168378e+00, -8.89181307e-17)|| = +6.33168378e+00 residual 2.40073360e-14
IRAM: Post Compression EigValue[0076]: ||(+6.46134095e+00, +7.36554440e-17)|| = +6.46134095e+00 residual 2.69335291e-14
IRAM: Post Compression EigValue[0077]: ||(+6.47683109e+00, -2.42370008e-17)|| = +6.47683109e+00 residual 3.01360898e-14
IRAM: Post Compression EigValue[0078]: ||(+6.70307632e+00, +5.85672461e-17)|| = +6.70307632e+00 residual 3.93748542e-14
IRAM: Post Compression EigValue[0079]: ||(+6.73923841e+00, +3.75828519e-18)|| = +6.73923841e+00 residual 2.44164094e-14
IRAM: Post Compression EigValue[0080]: ||(+6.89014735e+00, -1.95207213e-17)|| = +6.89014735e+00 residual 2.95075433e-14
IRAM: Post Compression EigValue[0081]: ||(+6.96401323e+00, +9.17985721e-17)|| = +6.96401323e+00 residual 3.74393525e-14
IRAM: Post Compression EigValue[0082]: ||(+7.20994999e+00, -7.02232665e-17)|| = +7.20994999e+00 residual 3.90727227e-14
IRAM: Post Compression EigValue[0083]: ||(+7.21831072e+00, -1.55891755e-16)|| = +7.21831072e+00 residual 3.19811426e-14
IRAM: Post Compression EigValue[0084]: ||(+7.41295610e+00, +1.22389485e-16)|| = +7.41295610e+00 residual 3.48496506e-14
IRAM: Post Compression EigValue[0085]: ||(+7.47203720e+00, -8.10187014e-17)|| = +7.47203720e+00 residual 3.35585534e-14
IRAM: Post Compression EigValue[0086]: ||(+7.68588751e+00, +4.63432901e-18)|| = +7.68588751e+00 residual 4.75208753e-14
IRAM: Post Compression EigValue[0087]: ||(+7.69280102e+00, +5.30259566e-17)|| = +7.69280102e+00 residual 3.30237090e-14
IRAM: Post Compression EigValue[0088]: ||(+7.97301663e+00, +2.21329709e-17)|| = +7.97301663e+00 residual 3.63596538e-14
IRAM: Post Compression EigValue[0089]: ||(+8.02465232e+00, -2.37313221e-17)|| = +8.02465232e+00 residual 3.44970393e-14
IRAM: Post Compression EigValue[0090]: ||(+8.24103952e+00, -3.69755292e-17)|| = +8.24103952e+00 residual 4.22563281e-14
IRAM: Post Compression EigValue[0091]: ||(+8.26625360e+00, -3.19331421e-18)|| = +8.26625360e+00 residual 3.29142112e-14
IRAM: Post Compression EigValue[0092]: ||(+8.58205939e+00, -7.08119544e-17)|| = +8.58205939e+00 residual 3.30281017e-14
IRAM: Post Compression EigValue[0093]: ||(+8.58383048e+00, -1.39328582e-17)|| = +8.58383048e+00 residual 2.80543623e-14
IRAM: Post Compression EigValue[0094]: ||(+8.77455568e+00, -3.99941429e-17)|| = +8.77455568e+00 residual 4.77694484e-14
IRAM: Post Compression EigValue[0095]: ||(+8.78283446e+00, +1.04882478e-16)|| = +8.78283446e+00 residual 3.77268650e-14
IRAM: Post Compression EigValue[0096]: ||(+8.94542528e+00, -1.68728963e-18)|| = +8.94542528e+00 residual 4.86593415e-14
IRAM: Post Compression EigValue[0097]: ||(+8.96797202e+00, -2.51873717e-17)|| = +8.96797202e+00 residual 4.41951428e-14
IRAM: Post Compression EigValue[0098]: ||(+9.46488406e+00, -1.27539445e-16)|| = +9.46488406e+00 residual 3.40393033e-14
IRAM: Post Compression EigValue[0099]: ||(+9.58532672e+00, +1.07315686e-16)|| = +9.58532672e+00 residual 3.71156276e-14
IRAM: Post Compression EigValue[0100]: ||(+9.85161258e+00, +5.93945855e-17)|| = +9.85161258e+00 residual 5.56854495e-14
IRAM: Post Compression EigValue[0101]: ||(+9.94591340e+00, -1.16102806e-16)|| = +9.94591340e+00 residual 4.29679563e-14
IRAM: Post Compression EigValue[0102]: ||(+1.04019207e+01, -1.42189303e-16)|| = +1.04019207e+01 residual 5.42922251e-14
IRAM: Post Compression EigValue[0103]: ||(+1.04305851e+01, -2.63723708e-18)|| = +1.04305851e+01 residual 7.61921300e-14
IRAM: Post Compression EigValue[0104]: ||(+1.10271537e+01, -2.19415415e-17)|| = +1.10271537e+01 residual 6.21841215e-14
IRAM: Post Compression EigValue[0105]: ||(+1.10950302e+01, -4.99063342e-17)|| = +1.10950302e+01 residual 8.99236901e-14
IRAM: Post Compression EigValue[0106]: ||(+1.13732913e+01, -3.20601971e-18)|| = +1.13732913e+01 residual 8.80655708e-14
IRAM: Post Compression EigValue[0107]: ||(+1.14375919e+01, -9.78609986e-17)|| = +1.14375919e+01 residual 6.63323733e-14
IRAM: Post Compression EigValue[0108]: ||(+1.16845873e+01, +1.20473496e-16)|| = +1.16845873e+01 residual 4.47522561e-14
IRAM: Post Compression EigValue[0109]: ||(+1.16894539e+01, +1.53318257e-17)|| = +1.16894539e+01 residual 6.13111215e-14
IRAM: Post Compression EigValue[0110]: ||(+1.21169663e+01, +1.01643954e-20)|| = +1.21169663e+01 residual 6.76380912e-14
IRAM: Post Compression EigValue[0111]: ||(+1.21698762e+01, +2.72078841e-16)|| = +1.21698762e+01 residual 8.53591552e-14
IRAM: Post Compression EigValue[0112]: ||(+1.24142685e+01, -1.67837461e-17)|| = +1.24142685e+01 residual 8.37039126e-14
IRAM: Post Compression EigValue[0113]: ||(+1.24906008e+01, -5.42577542e-17)|| = +1.24906008e+01 residual 4.45362104e-14
IRAM: Post Compression EigValue[0114]: ||(+1.28211693e+01, -5.08316118e-17)|| = +1.28211693e+01 residual 6.92939634e-14
IRAM: Post Compression EigValue[0115]: ||(+1.28789109e+01, -7.77110377e-17)|| = +1.28789109e+01 residual 8.86060394e-14
IRAM: Post Compression EigValue[0116]: ||(+1.30176504e+01, -2.79117685e-16)|| = +1.30176504e+01 residual 9.36389161e-14
IRAM: Post Compression EigValue[0117]: ||(+1.30594773e+01, -3.53996244e-17)|| = +1.30594773e+01 residual 1.13534981e-13
IRAM: Post Compression EigValue[0118]: ||(+1.32423153e+01, -3.26073803e-17)|| = +1.32423153e+01 residual 5.85916716e-14
IRAM: Post Compression EigValue[0119]: ||(+1.32507605e+01, -8.09051990e-17)|| = +1.32507605e+01 residual 6.98551953e-14
IRAM: Post Compression EigValue[0120]: ||(+1.47553315e+01, -1.02429484e-16)|| = +1.47553315e+01 residual 9.44775862e-14
IRAM: Post Compression EigValue[0121]: ||(+1.47720854e+01, -9.14117957e-17)|| = +1.47720854e+01 residual 7.95005810e-14
IRAM: Post Compression EigValue[0122]: ||(+1.50895773e+01, -1.85364690e-16)|| = +1.50895773e+01 residual 1.08816418e-13
IRAM: Post Compression EigValue[0123]: ||(+1.51572583e+01, -2.04312817e-16)|| = +1.51572583e+01 residual 8.57592366e-14
IRAM: Post Compression EigValue[0124]: ||(+1.54360435e+01, -1.76785517e-16)|| = +1.54360435e+01 residual 1.15875996e-13
IRAM: Post Compression EigValue[0125]: ||(+1.54604865e+01, +1.77672360e-16)|| = +1.54604865e+01 residual 1.05692623e-13
IRAM: Post Compression EigValue[0126]: ||(+1.57060762e+01, -1.26894006e-17)|| = +1.57060762e+01 residual 1.77540436e-13
IRAM: Post Compression EigValue[0127]: ||(+1.57213820e+01, -3.56725385e-16)|| = +1.57213820e+01 residual 9.24497601e-14
IRAM: Algorithmic compression report: 
IRAM: Complex(double) elems pre = 16384 Complex(double) elems post = 17408
IRAM: Ratio: 106.25% of original data.
IRAM: 8 low eigenvectors used 
IRAM: 120 high eigenvectors reconstructed 
IRAM: Compress/decompress time = 0.012033
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 1, res = 5.3500382520434121e-13, truersq = 2.0273302319796679e-14
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 71, res = 1.7462162297479900e-08, truersq = 6.6097081440676235e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 71, res = 1.7337922204225600e-08, truersq = 6.5626811686016999e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 72, res = 2.0664003858301122e-08, truersq = 7.8216554566386063e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 72, res = 2.4389626880946524e-08, truersq = 9.2318640622615888e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 72, res = 2.3383166479876691e-08, truersq = 8.8509031452632758e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 72, res = 1.9632374660782942e-08, truersq = 7.4311694261426941e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 73, res = 1.6374788845549241e-08, truersq = 6.1981194362290985e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 73, res = 1.5246479352388583e-08, truersq = 5.7710366469429885e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 73, res = 1.5850459073894768e-08, truersq = 5.9996528527415861e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 73, res = 1.0893880966957245e-08, truersq = 4.1235086932785038e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 73, res = 1.0750075529217547e-08, truersq = 4.0690752434833479e-10
CG: Applying deflation space preconditioner
IRAM: Deflating 128 eigenpairs.
CG: Converged iter = 72, res = 2.4138940657169294e-08, truersq = 9.1369752348038132e-10
```

Here we see the output from one HMC integegration. The inverter calls the IRAM eigensolver
at the start of the HMC integration and computes a deflation space. It then uses this
deflation space for all CG calls throughout the integration. The user will notice
that the deflation space is effective for one CG call only, then its efficacy quickly
fades. This is because as the gauge field evolves, the dirac operator also evolves.

One of the goals of this code is to devise methods which will `evolve` the deflation
space with the gauge evolution so that the efficacy of the deflation is restored
throughout the HMC integration.

Notice also that the MG (compression) routines are also tested. These routines use low mode
coherence to reconstruct high lying eigenmodes from low eigenmodes.

## HMC evolution with spectrum inspection

In order to devise methods for eigenvector evolution, we need to observe how the eigenspace is
evolving during the HMC evolution. To do this, we set `--eig-inspection true`. Note, this is a
FLOP intense operation and will recompute the eigenspectrum each time the gauge field evolves.
and the inverter is called.

# Data files

The executable will create a variety of data files for analysis. We show here the eigensolver
related data files only from the previos run with `--eig-inspection true` and `--eig-deflate true`.

## Eigenvalue evolution

Here is how to inspect the evolution of the 0th eigenvalue:
```
grep " 0 " data/eig/eigenvalues*
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:0 0 1.4724738299611465e-01 -1.2709729373550777e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:1 0 1.3724051903580350e-01 1.3484764520288461e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:2 0 1.3767889341887238e-01 4.6677868140816356e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:3 0 1.2419298384207285e-01 9.9706365803671829e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:4 0 1.1359992377517858e-01 6.1760877953455275e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:5 0 1.1407900600646727e-01 5.5195578518846397e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:6 0 1.0975354069792145e-01 2.1300761041077518e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:7 0 1.0952649651778687e-01 -1.4608936059972528e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:8 0 1.0964447719413478e-01 1.4330738676358694e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:9 0 1.1343134676351144e-01 -1.8677076486957322e-18
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:10 0 1.2418375554653939e-01 1.8756485825762413e-19
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:11 0 1.2456016734026923e-01 -2.5318344189023071e-19
data/eig/eigenvalues100_LX8_LY8_B4.000000_MLight0.050000.dat:12 0 1.4532746673539923e-01 4.5900186016118501e-18
data/eig/eigenvalues101_LX8_LY8_B4.000000_MLight0.050000.dat:0 0 1.4532746673539901e-01 -7.2611899403374897e-18
data/eig/eigenvalues101_LX8_LY8_B4.000000_MLight0.050000.dat:1 0 1.4315436810112561e-01 -1.4442176448481838e-18
data/eig/eigenvalues101_LX8_LY8_B4.000000_MLight0.050000.dat:2 0 1.4287418400543220e-01 -5.8867731042989806e-18
data/eig/eigenvalues101_LX8_LY8_B4.000000_MLight0.050000.dat:3 0 1.3751999853066446e-01 -3.0456128409712437e-19
data/eig/eigenvalues101_LX8_LY8_B4.000000_MLight0.050000.dat:4 0 1.3240743046747494e-01 7.4186851289675861e-18
```

The first column is the call to the CG solver for that integration, the second colum is the eigenvalue index (0 being the smallest)
and the two floats are real and imaginary components respectively. One can inspect the Nth eigenvalue evolution
by using `grep " N " <filename>`

## Eigenvector evolution

Here is how to inspect the evolution of eigenvectors:
```
grep " 0 0 0 0 " data/eig/eigenvectors10*
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:0   0 0 0 0 -1.2449146401686168e-02 1.1727778252333600e-01 1.1793667588889178e-01 1.6765512333074728e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:1   0 0 0 0 2.7021671490851416e-02 -1.2171959838557347e-01 1.2468292329467093e-01 -1.3523397768840424e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:2   0 0 0 0 1.6279811368230316e-02 -1.2337206880924897e-01 1.2444155102078737e-01 -1.4395973039538541e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:3   0 0 0 0 -3.7846977769885296e-03 1.2882709152212518e-01 1.2888267318500637e-01 1.6001659995174424e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:4   0 0 0 0 -1.1091487453822985e-02 1.3452989773269122e-01 1.3498634922801148e-01 1.6530565425419488e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:5   0 0 0 0 -7.4578110571878881e-03 -1.3435640652725991e-01 1.3456323019563346e-01 -1.6262470935598055e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:6   0 0 0 0 -4.5081948550803289e-02 -1.3091050753447114e-01 1.3845556351432795e-01 -1.9024488068106471e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:7   0 0 0 0 -3.0828129617769293e-02 1.3667222198498929e-01 1.4010592363652644e-01 1.7926461109673026e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:8   0 0 0 0 7.1910328864788681e-02 1.1984030037553711e-01 1.3975976885907124e-01 1.0303391016848289e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:9   0 0 0 0 4.2463080988933408e-02 1.3437655899861059e-01 1.4092612552461586e-01 1.2647254334397644e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:10   0 0 0 0 1.4557573802524077e-02 -1.4280944730400849e-01 1.4354950781626635e-01 -1.4692101728958871e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:11   0 0 0 0 -3.8134005609884274e-02 1.3801374169465583e-01 1.4318517828467398e-01 1.8403761796899565e+00
data/eig/eigenvectors100_LX8_LY8_B4.000000_MLight0.050000.dat:12   0 0 0 0 4.5918071484391934e-02 -1.3593726901978009e-01 1.4348313628228856e-01 -1.2450413287691529e+00
data/eig/eigenvectors101_LX8_LY8_B4.000000_MLight0.050000.dat:0   0 0 0 0 -1.2798702659742443e-02 1.4291117383756335e-01 1.4348313628229192e-01 1.6601150893629613e+00
data/eig/eigenvectors101_LX8_LY8_B4.000000_MLight0.050000.dat:1   0 0 0 0 1.5906685507443244e-02 -1.2759282658684951e-01 1.2858052745363333e-01 -1.4467686777793300e+00
data/eig/eigenvectors101_LX8_LY8_B4.000000_MLight0.050000.dat:2   0 0 0 0 1.0487872349044613e-02 1.2874564022276139e-01 1.2917211519046426e-01 1.4895138467494771e+00
data/eig/eigenvectors101_LX8_LY8_B4.000000_MLight0.050000.dat:3   0 0 0 0 3.0354243571594722e-02 -1.1513755549017675e-01 1.1907156162181308e-01 -1.3130268320817393e+00
data/eig/eigenvectors101_LX8_LY8_B4.000000_MLight0.050000.dat:4   0 0 0 0 3.3279846724164775e-03 -1.1752241291413489e-01 1.1756952419373082e-01 -1.5424860203156072e+00
```

The first column is again the call to the CG solver. The next four columns are respectively
(eigenvector index) (x coord) (y coord) (spin coord). The two floats that follow are the
complex value of that component in (real,imag) cartesian, and the next two are the polar
(modulus, argument) values.

