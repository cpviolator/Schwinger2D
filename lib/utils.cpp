#include "utils.h"
#include "blas.h"

void Param::usage(char **argv) {

  printf("\n\n");
  printf("This is an exhaustive list of run-time options. Those not set at\n");
  printf("run time will be given a default value. Important input options will\n");
  printf("dumped to stdout at execution. Please ensure that they are sensible!\n\n");
  printf("Default values are given in parentheses\n\n");
  printf("KEY: <N>     = integer\n");
  printf("     <float> = decimal\n");
  printf("     <bool>  = true/false\n");
  
  // The user is advised to study these parameters and set them accordingly.
  //-----------------------------------------------------------------------
  printf("GENERAL PARAMS\n");
  printf("--help                           Print this message.\n");
  printf("--seed <N>                       Sets the RNG seed (1234).\n");
  printf("--verbosity <bool>               Sets the verbosity as either verbose or quiet (true).\n");
  printf("\nPHYSICS PARAMS\n");
  printf("--beta <float>                   Sets the beta value (3.0).\n");
  printf("--dim <N> <N>                    Sets the lattice dimensions (32 32).\n");
  printf("--mass <float>                   Sets the light, degenerate quark mass (0.1).\n");
  printf("--mass-heavy <float>             Sets the heavy quark mass (0.5).\n");
  printf("--flavours <N>                   Sets the number of fermion flavours in the simulation (2).\n");
  printf("                                 0=pure gauge, 2=two degenerate, 3=two degenerate + one heavy.\n");
  printf("\nHMC PARAMS\n");
  printf("--hmc-traj-length <float>        Sets the HMC trajectory length (1.0).\n");
  printf("--hmc-n-step <N>                 Sets the number of HMC outer steps (4).\n");
  printf("--hmc-inner-step <N>             Sets the number of HMC inner steps (1).\n");
  printf("--hmc-n-trajectories <N>         Sets the number of HMC trajectories after thermalisation (1000).\n");
  printf("--hmc-therm <N>                  Sets the number of HMC thermalisation steps (250).\n");
  printf("--hmc-checkpoint <N>             Checkpoint the gauge field with this frequency (100).\n");
  printf("--hmc-checkpoint-start <N>       Start the HMC from this saved checkpoint (0).\n");  
  printf("--hmc-skip <N>                   Skip this number of trajectories between measurements (5).\n");
  printf("--hmc-reverse <N>                Perform a reversibility check with this frequency (100).\n");  
  printf("--hmc-integrator <LEAPFROG/FGI>  Sets the HMC integrator (FGI).\n");
  printf("--pfe-degree <N>                 Degree of the rational polynomial required for heavy mass fermion (15).\n");
  printf("--pfe-prec <N>                   GMP bit-wise precision of the rational polynomal computation (50).\n");
  printf("--cg-max-iter <N>                Maximum CG iterations before failure exit (1000).\n");
  printf("--cg-tol <float>                 Relative residual norm of CG solution (1e-9).\n");
  printf("--cg-verbosity <bool>            Sets CG verbosity as verbose or quiet (false).\n");
  printf("\nEIGENSOLVER PARAMS\n");
  printf("--eig-n-ev <N>                   Size of IRAM search space (16).\n");
  printf("--eig-n-kr <N>                   Size of IRAM Krylov space (32).\n");
  printf("--eig-n-conv <N>                 Number of required converged eigenpairs (16).\n");
  printf("--eig-n-deflate <N>              Number of eigenpairs to use in CG deflation. (0)\n");
  printf("--eig-max-restarts <N>           Maximum number of IRAM restarts (100).\n");
  printf("--eig-tol <float>                Residual norm of eigenvectors (1e-9).\n");
  printf("--eig-operator <M,Mdag,MdagM, MMdag> Operator to eigensolve (MdagM).\n");
  printf("--eig-spectrum <SR, LR, SI, LI, SM, LM> Spectrum to compute (SR)\n"
	 "                                 (S)malest/(L)argest (R)eal/(I)maginary/(M)odulus.\n");
  printf("--eig-block-scheme <N> <N>       Size of coarsening scheme for MG projector (2 2).\n");
  printf("--eig-low-modes <N>              Number of low eigenmodes used to construct MG projector (16).\n");
  printf("--eig-verbosity <bool>           Sets IRAM verbosity as verbose or quiet (false).\n");
  printf("--eig-deflate <bool>             Compute a deflation space at the start of the HMC trajectory\n"
	 "                                 and use it throughout the HMC integration (false)\n");
  printf("--eig-feast <bool>               Use FEAST eigensolver (false: use Lanczos)\n");
  printf("--eig-feast-M0 <N>               M0 (search space size) value for FEAST (32) \n");
  printf("--eig-feast-Emax <float>         Maximum eigenvalue size for FEAST (1.0) \n");
  printf("--eig-feast-Ncontour <N>         Contour integration points for FEAST (8)\n");
  printf("--eig-feast-init-guess <bool>    Use previous eigenspace ats init guess in FEAST (true)\n");
  printf("--eig-inspection <bool>          Inspect the eigenspectrum at each call of CG (false).\n");
  printf("--eig-use-comp-space <bool>      Use the compressed space for deflation (false).\n");
  printf("\nMEASUREMENT PARAMS\n");
  printf("--smear-type <APE, WILSON>       Smearing algorithm (WILSON)\n");
  printf("--ape-alpha <float>              Projection coefficient for APE smearing (0.5).\n");
  printf("--ape-iter <N>                   Number of APE smearing hits (1).\n");
  printf("--wilson-time <float>            Wilson Flow time (1.0)\n");
  printf("--wilson-steps <N>               Wilson Flow steps (10)\n");  
  printf("--meas-pl <bool>                 Measure Polyakov loops every measurement interval (false).\n");
  printf("--meas-wl <bool>                 Measure Wilson loops every measurement interval (false).\n");
  printf("--meas-pc <bool>                 Measure Pion every measurement interval (false).\n");
  printf("--meas-vt <bool>                 Measure Vacuum trace every measurement interval (false).\n");
}

int Param::init(int argc, char **argv, int *idx) {

  int ret = -1;
  int i = *idx;

  // Get help!
  if( strcmp(argv[i], "--help") == 0){
    usage(argv);
    exit(0);
  }

  // RNG seed
  if( strcmp(argv[i], "--seed") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    seed = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Verbosity
  if( strcmp(argv[i], "--verbosity") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string verbose(argv[i+1]);
    if (verbose == "V" || verbose == "v" ||
	verbose == "yes" || verbose == "YES" ||
	verbose == "true" || verbose == "TRUE" || verbose == "1" ) {
      verbosity = true;
      eig_param.verbosity = true;
    }
    else if (verbose == "Q" || verbose == "q" ||
	     verbose == "no" || verbose == "NO" ||
	     verbose == "false" || verbose == "FALSE" || verbose == "0") {
      verbosity = false;
      eig_param.verbosity = false;
    } else {
      cout<<"Invalid verbosity condition ("<< verbose << ") given. Use v/q for verbose/quiet."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // Physics parameters
  //-------------------------
  // Beta coupling
  if( strcmp(argv[i], "--beta") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    beta = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Lattice dim size
  if( strcmp(argv[i], "--dim") == 0){
    if (i+2 >= argc){
      usage(argv);
    }
    Nx = atoi(argv[i+1]);
    Ny = atoi(argv[i+2]);
    eig_param.Nx = Nx;
    eig_param.Ny = Ny;
    i+=2;
    ret = 0;
    goto out;
  }

  // Light degenerate mass
  if( strcmp(argv[i], "--mass") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    m = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Heavy mass
  if( strcmp(argv[i], "--mass-heavy") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    m_heavy = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Fermion flavours
  if( strcmp(argv[i], "--flavours") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    flavours = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // HMC params
  //-------------------------
  // HMC trajectory length
  if( strcmp(argv[i], "--hmc-traj-length") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    tau = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Outer HMC steps
  if( strcmp(argv[i], "--hmc-n-step") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_step = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Inner HMC steps
  if( strcmp(argv[i], "--hmc-inner-step") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    inner_step = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Number HMC of thermalised HMC trajectories
  if( strcmp(argv[i], "--hmc-n-trajectories") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    iter_hmc = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Number HMC of thermalisation steps
  if( strcmp(argv[i], "--hmc-therm") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    therm = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Checkpoint 
  if( strcmp(argv[i], "--hmc-checkpoint") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    chkpt = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Checkpoint start
  if( strcmp(argv[i], "--hmc-checkpoint-start") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    checkpoint_start = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // HMC skip
  if( strcmp(argv[i], "--hmc-skip") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    skip = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Reversibility check frequency
  if( strcmp(argv[i], "--hmc-reverse") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    reverse = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // HMC Integrator
  if( strcmp(argv[i], "--hmc-integrator") == 0){
    if (i+1 >= argc){
      usage(argv);
    }    
    std::string hmc_integrator(argv[i+1]);
    if (hmc_integrator == "LEAPFROG" || hmc_integrator == "leapfrog" || hmc_integrator == "0") {
      integrator = LEAPFROG;
    } else if(hmc_integrator == "FGI" || hmc_integrator == "fgi" || hmc_integrator == "1") {
      integrator  = FGI;
    } else {
      cout<<"Invalid integrator ("<< hmc_integrator <<") given. Use LEAPFROG or FGI"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Degree of partial fraction expansion
  if( strcmp(argv[i], "--pfe-degree") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    pfe_degree = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Bit-wise precision of partial fraction expansion
  if( strcmp(argv[i], "--pfe-prec") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    pfe_prec = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Maximun CG iteration  
  if( strcmp(argv[i], "--cg-max-iter") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    max_iter_cg = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Tolerance on CG residual norm
  if( strcmp(argv[i], "--cg-tol") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    tol_cg = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // CG verbosity
  if( strcmp(argv[i], "--cg-verbosity") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string cg_verbose(argv[i+1]);
    if (cg_verbose == "V" || cg_verbose == "v" ||
	cg_verbose == "yes" || cg_verbose == "YES" ||
	cg_verbose == "true" || cg_verbose == "TRUE" ||
	cg_verbose == "1") {
      cg_verbosity = true;
    }
    else if (cg_verbose == "Q" || cg_verbose == "q" ||
	     cg_verbose == "no" || cg_verbose == "NO" ||
	     cg_verbose == "false" || cg_verbose == "FALSE" ||
	     cg_verbose == "0") {
      cg_verbosity = false;
    } else {
      cout<<"Invalid CG verbosity condition ("<< cg_verbose<< ") given. Use v/q for verbose/quiet"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // Eigensolver params
  //-------------------------------------
  // Size of search space
  if( strcmp(argv[i], "--eig-n-ev") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.n_ev = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Size of Krylov space
  if( strcmp(argv[i], "--eig-n-kr") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.n_kr = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Number of converged vectors desired
  if( strcmp(argv[i], "--eig-n-conv") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.n_conv = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Number of vectors to deflate
  if( strcmp(argv[i], "--eig-n-deflate") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.n_deflate = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Maximum number of restarts
  if( strcmp(argv[i], "--eig-max-restarts") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.max_restarts = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Tolerance on eigenvector residual
  if( strcmp(argv[i], "--eig-tol") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.tol = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Spectrum to compute
  if( strcmp(argv[i], "--eig-spectrum") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string spec(argv[i+1]);
    if (spec == "SR" || spec == "sr") eig_param.spectrum = SR;
    else if (spec == "LR" || spec == "lr") eig_param.spectrum = LR;
    else if (spec == "SI" || spec == "si") eig_param.spectrum = SI;
    else if (spec == "LI" || spec == "li") eig_param.spectrum = LI;
    else if (spec == "SM" || spec == "sm") eig_param.spectrum = SM;
    else if (spec == "LM" || spec == "lm") eig_param.spectrum = LM;
    else {
      cout<<"Invalid spectrum ("<< spec <<") given. Use SR, LR, SI, LI, SM, LM" << endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Eigensolver verbosity
  if( strcmp(argv[i], "--eig-verbosity") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_verbose(argv[i+1]);
    if (eig_verbose == "V" || eig_verbose == "v" ||
	eig_verbose == "yes" || eig_verbose == "YES" ||
	eig_verbose == "true" || eig_verbose == "TRUE" ||
	eig_verbose == "1") {
      eig_param.iram_verbose = true;
    }
    else if (eig_verbose == "Q" || eig_verbose == "q" ||
	     eig_verbose == "no" || eig_verbose == "NO" ||
	     eig_verbose == "false" || eig_verbose == "FALSE" ||
	     eig_verbose == "0") {
      eig_param.iram_verbose = false;
    } else {
      cout<<"Invalid Eigensolver verbosity condition ("<< eig_verbose<< ") given. Use v/q for verbose/quiet"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Operator to eigensolve
  if( strcmp(argv[i], "--eig-operator") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_op(argv[i+1]);
    if (eig_op == "M" || eig_op == "m") eig_param.op = M;
    else if (eig_op == "Mdag" || eig_op == "mdag") eig_param.op = Mdag;
    else if (eig_op == "MdagM" || eig_op == "mdagm") eig_param.op = MdagM;
    else if (eig_op == "MMdag" || eig_op == "mmdag") eig_param.op = MMdag;
     else {
      cout<<"Invalid eigen operator ("<< eig_op <<") given. Use M, Mdag, MdagM, MMdag" << endl;
      exit(0);
     }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Eig MG block scheme
  if( strcmp(argv[i], "--eig-block-scheme") == 0){
    if (i+2 >= argc){
      usage(argv);
    }
    eig_param.block_scheme[0] = atoi(argv[i+1]);
    eig_param.block_scheme[1] = atoi(argv[i+2]);
    i+=2;
    ret = 0;
    goto out;
  }
  
  // Number of low modes to use in MG 
  if( strcmp(argv[i], "--eig-low-modes") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.n_low = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Number of contour integration points
  if( strcmp(argv[i], "--eig-feast-Ncontour") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.feast_Ncontour = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // M0 search space size
  if( strcmp(argv[i], "--eig-feast-M0") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.feast_M0 = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Tolerance on eigenvector residual
  if( strcmp(argv[i], "--eig-feast-Emax") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    eig_param.feast_Emax = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  // Eigensolver FEAST
  if( strcmp(argv[i], "--eig-feast") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_feast(argv[i+1]);
    if (eig_feast == "yes" || eig_feast == "YES" ||
	eig_feast == "true" || eig_feast == "TRUE" ||
	eig_feast == "1") {
      use_feast = true;
    }
    else if (eig_feast == "no" || eig_feast == "NO" ||
	     eig_feast == "false" || eig_feast == "FALSE" ||
	     eig_feast == "0") {
      use_feast = false;
    } else {
      cout<<"Invalid Eigensolver FEAST condition ("<< eig_feast << ") given. Use true/false"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  // Eigensolver FEAST init guess
  if( strcmp(argv[i], "--eig-feast-init-guess") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_feast_init(argv[i+1]);
    if (eig_feast_init == "yes" || eig_feast_init == "YES" ||
	eig_feast_init == "true" || eig_feast_init == "TRUE" ||
	eig_feast_init == "1") {
      eig_param.feast_init_guess = true;
    }
    else if (eig_feast_init == "no" || eig_feast_init == "NO" ||
	     eig_feast_init == "false" || eig_feast_init == "FALSE" ||
	     eig_feast_init == "0") {
      eig_param.feast_init_guess = false;
    } else {
      cout<<"Invalid Eigensolver FEAST init guess condition ("<< eig_feast_init << ") given. Use true/false"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Eigensolver deflation
  if( strcmp(argv[i], "--eig-deflate") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_defl(argv[i+1]);
    if (eig_defl == "yes" || eig_defl == "YES" ||
	eig_defl == "true" || eig_defl == "TRUE" ||
	eig_defl == "1") {
      deflate = true;
    }
    else if (eig_defl == "no" || eig_defl == "NO" ||
	     eig_defl == "false" || eig_defl == "FALSE" ||
	     eig_defl == "0") {
      deflate = false;
    } else {
      cout<<"Invalid Eigensolver deflation condition ("<< eig_defl << ") given. Use true/false"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // Eigensolver inspection
  if( strcmp(argv[i], "--eig-inspection") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_inspec(argv[i+1]);
    if (eig_inspec == "yes" || eig_inspec == "YES" ||
	eig_inspec == "true" || eig_inspec == "TRUE" ||
	eig_inspec == "1") {
      inspect_spectrum = true;
    }
    else if (eig_inspec == "no" || eig_inspec == "NO" ||
	     eig_inspec == "false" || eig_inspec == "FALSE" ||
	     eig_inspec == "0") {
      inspect_spectrum = false;
    } else {
      cout<<"Invalid Eigensolver inspection condition ("<< eig_inspec<< ") given. Use true/false"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // USe compressed space for deflation
  if( strcmp(argv[i], "--eig-use-comp-space") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string eig_ucs(argv[i+1]);
    if (eig_ucs == "yes" || eig_ucs == "YES" ||
	eig_ucs == "true" || eig_ucs == "TRUE" ||
	eig_ucs == "1") {
      eig_param.use_comp_space = true;
    }
    else if (eig_ucs == "no" || eig_ucs == "NO" ||
	     eig_ucs == "false" || eig_ucs == "FALSE" ||
	     eig_ucs == "0") {
      eig_param.use_comp_space = false;
    } else {
      cout<<"Invalid Use Compressed Space condition ("<< eig_ucs<< ") given. Use true/false"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  
  
  // Measurement params
  //-------------------------------------
  // HMC Integrator
  if( strcmp(argv[i], "--smear-type") == 0){
    if (i+1 >= argc){
      usage(argv);
    }    
    std::string smear(argv[i+1]);
    if (smear == "APE" || smear == "ape" || smear == "0") {
      smear_type = APE;
    } else if(smear == "WILSON" || smear == "wilson" || smear == "1") {
      smear  = WILSON;
    } else {
      cout<<"Invalid smear type ("<< smear <<") given. Use APE or WILSON"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 
  
  // APE smearing iterations 
  if( strcmp(argv[i], "--ape-iter") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    ape_iter = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // APE Alpha coefficient
  if( strcmp(argv[i], "--ape-alpha") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    ape_alpha = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Wilson Flow smearing iterations 
  if( strcmp(argv[i], "--wilson-steps") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    wilson_steps = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Wilson Flow time
  if( strcmp(argv[i], "--wilson-time") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    wilson_time = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  // Polyakov loops
  if( strcmp(argv[i], "--meas-pl") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string m_pl(argv[i+1]);
    if (m_pl == "yes" || m_pl == "YES" || m_pl == "true" || m_pl == "TRUE" || m_pl == "1") {
      meas_pl = true;
    }
    else if (m_pl == "no" || m_pl == "NO" || m_pl == "false" || m_pl == "FALSE" || m_pl == "0") {
      meas_pl = false;
    }
    else {
      cout<<"Invalid measure polyakov loop condition ("<< m_pl << ") given. Use true/false" << endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // Wilson loops 
  if( strcmp(argv[i], "--meas-wl") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string m_wl(argv[i+1]);
    if (m_wl == "yes" || m_wl == "YES" || m_wl == "true" || m_wl == "TRUE" || m_wl == "1") {
      meas_wl = true;
    }
    else if (m_wl == "no" || m_wl == "NO" || m_wl == "false" || m_wl == "FALSE" || m_wl == "0") {
      meas_wl = false;
    }
    else {
      cout<<"Invalid measure Wilson loop condition ("<< m_wl << ") given. Use true/false" << endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // Pseudo scalar correlation
  if( strcmp(argv[i], "--meas-pc") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string m_pc(argv[i+1]);
    if (m_pc == "yes" || m_pc == "YES" || m_pc == "true" || m_pc == "TRUE" || m_pc == "1") {
      meas_pc = true;
    }
    else if (m_pc == "no" || m_pc == "NO" || m_pc == "false" || m_pc == "FALSE" || m_pc == "0") {
      meas_pc = false;
    }
    else {
      cout<<"Invalid measure Pion condition ("<< m_pc << ") given. Use true/false" << endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  // Vacuum trace
  if( strcmp(argv[i], "--meas-vt") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string m_vt(argv[i+1]);
    if (m_vt == "yes" || m_vt == "YES" || m_vt == "true" || m_vt == "TRUE" || m_vt == "1") {
      meas_vt = true;
    }
    else if (m_vt == "no" || m_vt == "NO" || m_vt == "false" || m_vt == "FALSE" || m_vt == "0") {
      meas_vt = false;
    }
    else {
      cout<<"Invalid measure Vacuum trace condition ("<< m_vt << ") given. Use true/false" << endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  
  
 out:
  *idx = i;
  return ret ;
  
}

void Param::print() {
  
  cout<<endl;
  cout<<"********************************"<<endl;
  cout<<"*      Parameter status        *"<<endl;
  cout<<"********************************"<<endl;
  cout << endl;
  cout << "Physics:      X dim = "<< Nx << endl;
  cout << "              Y dim = "<< Ny << endl;
  cout << "              Beta = "<< beta << endl;
  cout << "              Flavours = " << flavours << endl;
  if (flavours > 0)
    cout << "              Mass = " << m << endl;
  if (flavours == 3)
    cout << "              Mass Heavy = " << m_heavy << endl;
  cout << "HMC:          Therm Sweeps = " << therm << endl; 
  cout << "              Data Points = " << iter_hmc << endl;
  cout << "              Start Point = " << checkpoint_start << endl;
  cout << "              Integrator = " << (integrator == LEAPFROG ? "LEAPFROG" : "FGI") << endl;
  cout << "              Trajectory Length = " << tau << endl;
  cout << "              Trajectory Steps = " << n_step << endl;
  if  (integrator == FGI) {
    cout << "              Inner Trajectory Steps = " << inner_step << endl;
  }
  cout << "Measurements: Polyakov loops = " << (meas_pl ? "true" : "false") << endl;
  cout << "              Wilson loops = " << (meas_wl ? "true" : "false") << endl;
  cout << "              Pion = " << (meas_pc ? "true" : "false") << endl;
  cout << "              Vacuum trace = " << (meas_vt ? "true" : "false") << endl;
  if(smear_type == APE) {
    cout << "              APE iter = " << ape_iter << endl;
    cout << "              APE alpha = " << ape_alpha << endl;
  } else {
    cout << "              Wilson Flow steps = " << wilson_steps << endl;
    cout << "              Wilson Flow time = " << wilson_time << endl;
  }
  if(deflate && flavours > 0) {
    cout << "Deflation:    nkv = " << eig_param.n_kr << endl;
    cout << "              nev = " << eig_param.n_ev << endl;
    cout << "              ndefl = " << eig_param.n_deflate << endl;
    cout << "              tol = " << eig_param.tol << endl;
    cout << "              maxiter = " << eig_param.max_restarts << endl;
  } else {
    cout << "Deflation:    OFF " << endl;
  }
}

void constructName(string &name, Param p) {

  // Basics
  name += "_LX" + to_string(p.Nx) + "_LY" + to_string(p.Ny) + "_B" + to_string(p.beta);

  // Fermions
  if(p.flavours == 2) {
    name += "_MLight"+ to_string(p.m);
  } else if(p.flavours == 3) {
    name += "_MLight"+ to_string(p.m);
    name += "_MHeavy"+ to_string(p.m_heavy);
  }
}

void gaussStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
	gauge->write(x, y, mu, polar(1.0,drand48()));
      }
}  

void coldStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, Complex(1.0,0.0));    
}  

// Normalized gaussian exp(-phi*phi/2) and <phi|phi> = 1
void gaussReal(field<double> *field) {
  
  double r, theta, sum;
  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      //for(int mu=0; mu<2; mu++) {	
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field->write(x,y,0, r*cos(theta));
      field->write(x,y,1, r*sin(theta));
    }
  }
}  


//normalized gaussian exp(-eta*eta/2) and <eta|eta> = 1;
void gaussComplex(field<Complex> *field) {

  double r1, theta1, r2, theta2, sum;
  double inv_sqrt2 = 1.0/sqrt(2);

  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	
	field->write(x,y,mu, Complex(r1*cos(theta1),r1*sin(theta1))*inv_sqrt2);
      }
    }
  }
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(field<Complex> *smeared, field<Complex> *gauge){

  blas::copy(smeared, gauge);

  if(gauge->p.smear_type == WILSON) {
    double wilson_flow_tau = gauge->p.wilson_time;
    int wilson_flow_steps = gauge->p.wilson_steps;
    double dt = wilson_flow_tau/wilson_flow_steps;
    for(int i=0; i<wilson_flow_steps; i++) wilsonFlow(smeared, dt);
  } else if(gauge->p.smear_type == APE) {

    double alpha = gauge->p.ape_alpha;
    Complex tmp;
    
    field<Complex> *smeared_tmp = new field<Complex>(gauge->p);
    
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
    for(int i=0; i<gauge->p.ape_iter; i++) {
      blas::copy(smeared_tmp, smeared);
      
      for(int x=0; x<Nx; x++) {
	int xp1 = (x+1)%Nx;
	int xm1 = (x-1+Nx)%Nx;
	for(int y=0; y<Ny; y++) {
	  int yp1 = (y+1)%Ny;
	  int ym1 = (y-1+Ny)%Ny;
	  
	  //|->-|   |   |
	  //^   v + v   ^
	  //|   |   |->-|
	  tmp = smeared->read(x,y,0);
	  tmp += alpha * (smeared->read(x,y,1) * smeared->read(x,yp1,0) * conj(smeared->read(xp1,y,1)));
	  tmp += alpha * (conj(smeared->read(x,ym1,1)) * smeared->read(x,ym1,0) * smeared->read(xp1,ym1,1));			
	  smeared_tmp->write(x,y,0, tmp);
	  
	  //|->-    -<-|
	  //^    +     ^
	  //|-<-    ->-|
	  tmp = smeared->read(x,y,1);
	  tmp += alpha * (smeared->read(x,y,0) * smeared->read(xp1,y,1) * conj(smeared->read(x,yp1,0)));
	  tmp += alpha * (conj(smeared->read(xm1,y,0)) * smeared->read(xm1,y,1) * smeared->read(xm1, yp1,01));
	  smeared_tmp->write(x,y,1, tmp);	
	}
      }
      
      //Project back to U(1)
      for(int x=0; x<Nx; x++)
	for(int y=0; y<Ny; y++)
	  for(int mu=0; mu<2; mu++)
	    smeared->write(x,y,mu, polar(1.0,arg(smeared_tmp->read(x,y,mu))));
    }
    
    delete smeared_tmp;
  } else {
    cout << "Error Smear: Unknown smear type passed " << endl;
    exit(0);
  }
}

// staple x is 0th, y is 1st.
void wilsonFlow(field<Complex> *gauge, double dt) {

  // perform one time step of Runga Kutta integration for the Wilson flow
  // algorithm.

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;

  field<double> *z0 = new field<double>(gauge->p);
  field<double> *z1 = new field<double>(gauge->p);
  field<double> *z2 = new field<double>(gauge->p);
  
  field<Complex> *w0 = new field<Complex>(gauge->p);
  field<Complex> *w1 = new field<Complex>(gauge->p);
  field<Complex> *w2 = new field<Complex>(gauge->p);
  field<Complex> *w3 = new field<Complex>(gauge->p);
  
  blas::copy(w0, gauge);

  for (int x = 0; x < Nx; x++) {
    for (int y = 0; y < Ny; y++) {
      for (int mu = 0; mu < 2; mu++) {
	Complex s = staple(w0, x, y, mu);
	double z = -imag(w0->read(x, y, mu) * conj(s));
	z0->write(x, y, mu, dt * z);
	w1->write(x, y, mu, polar(1.0, (1.0/4.0) * z0->read(x, y, mu)) * w0->read(x, y, mu));
      }
    }
  }

  for (int x = 0; x < Nx; x++) {
    for (int y = 0; y < Ny; y++) {
      for (int mu = 0; mu < 2; mu++) {
	Complex s = staple(w1, x, y, mu);
	double z = -imag(w1->read(x, y, mu) * conj(s));
	z1->write(x, y, mu, dt * z);
	w2->write(x, y, mu, polar(1.0, (8.0/9.0) * z1->read(x, y, mu) - (17.0/36.0) * z0->read(x, y, mu)) * w1->read(x, y, mu));
      }
    }
  }

  for (int x = 0; x < Nx; x++) {
    for (int y = 0; y < Ny; y++) {
      for (int mu = 0; mu < 2; mu++) {
	Complex s = staple(w2, x, y, mu);
	double z = -imag(w2->read(x, y, mu) * conj(s));
	z2->write(x, y, mu, dt * z);
	w3->write(x, y, mu, polar(1.0, (3.0/4.0) * z2->read(x, y, mu) - (8.0/9.0) * z1->read(x, y, mu) + (17.0/36.0) * z0->read(x, y, mu)) * w2->read(x, y, mu));
      }
    }
  }
  
  blas::copy(gauge, w3);
}

Complex staple(field<Complex> *gauge, int x, int y, int mu) {
  
  // calculate the sum of staples attached to the link at (x,y,mu) in the
  // reverse direction as compared to the direction of the link
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  int xp1 = (x + 1) % Nx;
  int xm1 = (x - 1 + Nx) % Nx;
  int yp1 = (y + 1) % Ny;
  int ym1 = (y - 1 + Ny) % Ny;

  if (mu == 0) {
    return gauge->read(x,y,1) * gauge->read(x,yp1,0) * conj(gauge->read(xp1,y,1))
      + conj(gauge->read(x,ym1,1)) * gauge->read(x,ym1,0) * gauge->read(xp1,ym1,1);
  } else {
    return gauge->read(x,y,0) * gauge->read(xp1,y,1) * conj(gauge->read(x,yp1,0))
      + conj(gauge->read(xm1,y,0)) * gauge->read(xm1,y,1) * gauge->read(xm1,yp1,0);
  }
}
