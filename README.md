# Build

Make a build directory and enter it. Then invoke the CMake CL, configure the build, and make. For a quick set up:

```
1. mkdir build
2. cmake /path/to/Schwinger2D/
3. make -j
```

There is a launcher script in the `wilson` directory. Place the launcher script and the executable `wilson2D`
in the same directory in which you wish to run and launch with `./launcher.sh`

# Example output

Here is some example output from Mac OSX:

```
./wilson2D 3.0 1000 10 5 50 0 3 1.0 1 0.5 1234 1 0.1 10000 1e-18 0 56 24 24 1e-10 1000 0 11 1.0 100 4 4 16 16 1 1 8 8 1 0.6 1 2 100 12 40 0
Warning: requested Wilson loop max 16 greater than 4, truncating.

Physics:  XSize = 8
          YSize = 8
          Beta = 3
          Dynamic = True
          Flavours = 2
          Mass = 0.1
HMC:      Therm Sweeps: (10 accept) (10 accept/reject)
          Data Points = 1000
          Integrator = FGI
          Time Step = 0.333333
          Trajectory Steps = 3
          Inner Trajectory Steps = 1
          Trajectory Length = 1
Smearing: APE iter = 1
          APE alpha = 0.5
1 0.0760410000000000 
2 0.1360320000000000 
3 0.1956350000000000 
4 0.2565900000000000 
5 0.3157250000000000 
6 0.3765570000000000 
7 0.4383700000000000 
8 0.5002180000000001 
9 0.5610080000000000 
10 0.6218610000000000 
11 0.6947590000000000 
12 0.7665419999999999 
13 0.8396100000000000 
14 0.9118510000000000 
15 0.9839510000000000 
16 1.0561570000000000 
17 1.1277390000000000 
18 1.1997560000000000 
19 1.2716499999999999 
20 1.3427650000000000 
25 1.7032130000000001 0.7921973532623939 0.7921973532623939 1.2290477088449516 0.9824941355922718 -0.2062396490696301 0.0329456760470123 1.0000000000000000 0
30 2.0776289999999999 0.8375559304531049 0.8148766418577494 1.0822655045701117 0.9179701513176296 -0.0790565334744429 -0.0066494900850415 0.9000000000000000 1
35 2.4518840000000002 0.8448225294628630 0.8248586043927872 1.6910992254279704 1.0215375755283560 -0.5253787467448774 -0.0629968894415867 0.9333333333333333 0
40 2.8300280000000000 0.8240904034860990 0.8246665541661151 1.0265588306352749 0.8625274644250196 -0.0262122677052901 -0.0433417184636994 0.8000000000000000 0
45 3.2123400000000002 0.8412433932397022 0.8279819219808326 0.3801492347914480 0.7636609593211664 +0.9671913802218626 -0.0273943930943767 0.7200000000000000 -1
...
```

Points to note:
1. The command is printed to the terminal
2. The physics and algorithmic parameters are printed.
3. 10 thermalisation sweeps are performed and always accepted.
4. 10 further thermalisation sweeps are performed with metropolis.
5. After the thermalisation sweeps are performed, N_SKIP sweeps are
   performed and then measurements are taken.
6. From left to right, the data in the terminal is
   (wall clock) (lattice action) (average lattice action) (exp(-dH)) (dH) (average dH) (Average acceptance) (topological charge)

This data is dumped to file `filename.dat` in a directory data/data/filename.dat, where `filename` is constructed
automatically, e.g., `data_LX8_LY8_B3.000000_M0.100000_tau1.000000_nHMCstep3.dat `

Other data collected are in the directory `data` such as topological charge and pseudoscalar correlation data
