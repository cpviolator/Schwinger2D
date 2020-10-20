#include "utils.h"

void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  XSize = "<< p.Nx << endl;
  cout << "          YSize = "<< p.Ny << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  if (p.dynamic == true) cout << "          Mass = " << p.m << endl;
  cout << "HMC:      Therm Sweeps: (" << p.therm << " accept) (" << p.therm << " accept/reject)" << endl; 
  cout << "          Data Points = " << p.iter_hmc << endl;
  cout << "          Time Step = " << p.tau/p.n_step << endl;
  cout << "          Trajectory Steps " << p.n_step << endl;
  cout << "          Trajectory Length = " << p.tau << endl;
  cout << "Smearing: APE iter = " << p.smear_iter << endl;
  cout << "          APE alpha = " << p.alpha << endl;
#ifdef USE_VOATOL
  cout << "VOATOL:   nkv = " << p.n_kr << endl;
  cout << "          nev = " << p.n_ev << endl;
  cout << "          tol = " << p.eig_tol << endl;
  cout << "          maxiter = " << p.eig_max_restart << endl;
#endif  
}

void constructName(string &name, param_t p) {
  name += "_LX" + to_string(p.Nx) + "_LY" + to_string(p.Ny) + "_B" + to_string(p.beta);
  if(p.dynamic == true) name += "_M"+ to_string(p.m);
  name += "_tau" + to_string(p.tau) + "_nHMCstep" + to_string(p.n_step);
}

void gaussStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, polar(1.0,drand48()));  
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
      for(int mu=0; mu<2; mu++) {
	
	r = sqrt(-2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field->write(x,y,mu, r*cos(theta));
      }
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
	r2 = sqrt(-2.0*log(drand48()));
	theta2 = TWO_PI*(drand48());
	
	field->write(x,y,mu, Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2);
      }
    }
  }
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(field<Complex> *smeared, field<Complex> *gauge){

  double alpha = gauge->p.alpha;
  int iter = gauge->p.smear_iter;
  Complex tmp = 0;
  field<Complex> *smeared_tmp = new field<Complex>(gauge->p);
  smeared->copy(gauge);
  smeared_tmp->copy(smeared);
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int i=0; i<iter; i++) {    
    for(int x=0; x<Nx; x++) {
      for(int y=0; y<Ny; y++) {

	int xp1 = (x+1)%Nx;
	int xm1 = (x-1+Nx)%Nx;
	int yp1 = (y+1)%Ny;
	int ym1 = (y-1+Ny)%Ny;
		
	//|->-|   |   |
	//^   v + v   ^
	//|   |   |->-|
	tmp = alpha * (smeared->read(x,y,1) * smeared->read(x,yp1,0) * conj(smeared->read(xp1,y,1)));
	
	tmp += alpha * (conj(smeared->read(x,ym1,1)) * smeared->read(x,ym1,0) * smeared->read(xp1,ym1,1));
			
	smeared_tmp->write(x,y,0, tmp);
				
	//|->-    -<-|
	//^    +     ^
	//|-<-    ->-|
	tmp = alpha * (smeared->read(x,y,0) * smeared->read(xp1,y,1) * conj(smeared->read(x,yp1,0)));
	
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
}

void writeBlockToVector(std::vector<field<Complex> *> &kSpace,
			std::vector<std::vector<Complex>> &block_data_ortho,
			std::vector<std::vector<Complex>> &block_coef,
			int blockScheme[2], int nLow, int nConv) {
    
}

// Read the block data from the (iEig)th vector in kSpace 
void readVectorToBlock(std::vector<field<Complex> *> &kSpace,
		       std::vector<std::vector<Complex>> &block_data,
		       int blockScheme[2], int iEig) {

  // This will become class data
  int Nx = kSpace[0]->p.Nx;
  int Ny = kSpace[0]->p.Ny;
  if(Nx%blockScheme[0] != 0) {
    cout << "Error: x blockScheme = " << blockScheme[0] << " does not divide Nx = " << Nx << endl;
  }
  if(Ny%blockScheme[1] != 0) {
    cout << "Error: y blockScheme = " << blockScheme[1] << " does not divide Ny = " << Ny << endl;
  }
  
  int x_block_size = Nx/blockScheme[0];
  int y_block_size = Ny/blockScheme[1];
  int n_blocks = blockScheme[0]*blockScheme[1];
  int blk_size = 2 * x_block_size * y_block_size;
  
  for(int by=0; by<blockScheme[1]; by++) {
    for(int bx=0; bx<blockScheme[0]; bx++) {
      int blk_idx = by * blockScheme[0] + bx;
      
      // Location of the start of the desired block (x runs fastest)
      int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));
      
      for(int nx=0; nx<x_block_size; nx++) {
	for(int ny=0; ny<y_block_size; ny++) {
	  for(int mu=0; mu<2; mu++) {
	    
	    int loc_idx = 2*(nx + x_block_size * ny) + mu;      // Local index
	    int glo_idx = blk_offset + 2*((ny * Nx) + nx) + mu; // Global index
	    block_data[blk_idx][blk_size * iEig + loc_idx] = kSpace[iEig]->data[glo_idx];
	  }
	}
      }
    }
  }
}



void blockCompress(std::vector<field<Complex> *> &kSpace,
		   std::vector<std::vector<Complex>> &block_data_ortho,
		   std::vector<std::vector<Complex>> &block_coef,
		   int blockScheme[2], int nLow, int nConv) {
  
  int Nx = kSpace[0]->p.Nx;
  int Ny = kSpace[0]->p.Ny;
  if(Nx%blockScheme[0] != 0) {
    cout << "Error: x blockScheme = " << blockScheme[0] << " does not divide Nx = " << Nx << endl;
  }
  if(Ny%blockScheme[1] != 0) {
    cout << "Error: y blockScheme = " << blockScheme[1] << " does not divide Ny = " << Ny << endl;
  }
  
  int x_block_size = Nx/blockScheme[0];
  int y_block_size = Ny/blockScheme[1];
  int n_blocks = blockScheme[0]*blockScheme[1];
  int blk_size = 2 * x_block_size * y_block_size;

  // Object to hold the blocked eigenvector data
  std::vector<std::vector<Complex>> block_data(n_blocks, std::vector<Complex> (nConv * blk_size, 0.0));
  
  // Copy data from the eigenvector array into a block array
  //#pragma omp parallel for collapse(2)
  
  for(int i=0; i<nConv; i++) {
    readVectorToBlock(kSpace, block_data, blockScheme, i);
  }
  
  /*  
  for(int i=0; i<nConv; i++) {
    for(int by=0; by<blockScheme[1]; by++) {
      for(int bx=0; bx<blockScheme[0]; bx++) {
	int blk_idx = by * blockScheme[0] + bx;
	
	// Location of the start of the desired block (x runs fastest)
	int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));
	
	for(int nx=0; nx<x_block_size; nx++) {
	  for(int ny=0; ny<y_block_size; ny++) {
	    for(int mu=0; mu<2; mu++) {
	      
	      int loc_idx = 2*(nx + x_block_size * ny) + mu;      // Local index
	      int glo_idx = blk_offset + 2*((ny * Nx) + nx) + mu; // Global index
	      block_data[blk_idx][blk_size * i + loc_idx] = kSpace[i]->data[glo_idx];
	    }
	  }
	}
      }
    }
  }
  */
  
  // Compute an orthonormal basis from the nLow vectors.
#pragma omp parallel for collapse(2)
  for(int by=0; by<blockScheme[1]; by++) {
    for(int bx=0; bx<blockScheme[0]; bx++) {
      int blk_idx = by * blockScheme[0] + bx;
      
      // Location of the start of the desired block (x runs fastest)
      int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));
      
      for(int i=0; i<nLow; i++) {
	
	// Copy data from the block data into the block ortho array
	for(int k=0; k<blk_size; k++) {
	  block_data_ortho[blk_idx][blk_size * i + k] = block_data[blk_idx][blk_size * i + k];
	}
	
	// Loop up to i
	for(int j=0; j<i; j++) {
	  
	  // Copy data from the jth ortho array into a temp array
	  std::vector<Complex> temp(blk_size, 0.0);
	  for(int k=0; k<blk_size; k++) temp[k] = block_data_ortho[blk_idx][blk_size * j + k];
	  
	  // <Vj|Vi> : inner product
	  Complex ip = 0.0;	  
	  for(int k=0; k<blk_size; k++) ip += conj(temp[k]) * block_data[blk_idx][blk_size * i + k];
	  
	  // |Vi> = |Vi> - <Vj|Vi>|Vj> : CAXPY project and write to block_data_ortho 
	  for(int k=0; k<blk_size; k++) {
	    block_data_ortho[blk_idx][blk_size * i + k] = (block_data_ortho[blk_idx][blk_size * i + k] - ip * temp[k]);
	  }	  
	}

	// Normalize
	Complex ip = 0.0;
	for(int k=0; k<blk_size; k++) {
	  ip += conj(block_data_ortho[blk_idx][blk_size * i + k]) * block_data_ortho[blk_idx][blk_size * i + k];
	}
	
	double nrm = 1.0/sqrt(ip.real());
	for(int k=0; k<blk_size; k++) {
	  block_data_ortho[blk_idx][blk_size * i + k] *= nrm;
	}
      }
    }
  }
  
  // Get coefficients: project the blocked nConv eigenvectors on to the
  // orthonormalised low blocks.
  // Modified Gramm-Schmidt
#pragma omp parallel for collapse(3)
  for(int j=0; j<nConv; j++) {
    for(int by=0; by<blockScheme[1]; by++) {
      for(int bx=0; bx<blockScheme[0]; bx++) {
	
	int blk_idx = by * blockScheme[0] + bx;        
	for(int i=0; i<nLow; i++) {
	  

	  // Inner product between orthed i block (low) and j block (high)
	  Complex ip = 0.0;
	  for(int k=0; k<blk_size; k++) {
	    ip += (conj(block_data_ortho[blk_idx][blk_size * i + k]) * 
		   block_data[blk_idx][blk_size * j + k]);
	  }
	  
	  block_coef[blk_idx][nLow * j + i] = ip;
	  
	  for(int k=0; k<blk_size; k++) {
	    block_data[blk_idx][blk_size * j + k] = block_data[blk_idx][blk_size * j + k] - (ip * block_data_ortho[blk_idx][blk_size * i + k]);
	  }
	}
      }
    }
  }  
}

void blockExpand(std::vector<field<Complex> *> &kSpace,
		 std::vector<std::vector<Complex>> &block_data_ortho,
		 std::vector<std::vector<Complex>> &block_coef,
		 int blockScheme[2], int nLow, int nConv) {

  int Nx = kSpace[0]->p.Nx;
  int Ny = kSpace[0]->p.Ny;
  if(Nx%blockScheme[0] != 0) {
    cout << "Error: x blockScheme = " << blockScheme[0] << " does not divide Nx = " << Nx << endl;
  }
  if(Ny%blockScheme[1] != 0) {
    cout << "Error: y blockScheme = " << blockScheme[1] << " does not divide Ny = " << Ny << endl;
  }
  int x_block_size = Nx/blockScheme[0];
  int y_block_size = Ny/blockScheme[1];
  int n_blocks = blockScheme[0]*blockScheme[1];
  int blk_size = 2 * x_block_size * y_block_size;

  // Loop over the desired eigenvalues.
#pragma omp parallel for 
  for(int j=0; j<nConv; j++) {    
    std::vector<std::vector<Complex>> block_data_temp(n_blocks, std::vector<Complex> (blk_size, 0.0));      
    // Loop over blocks, Gramm-Schmidt the blocks on the low modes
    for(int by=0; by<blockScheme[1]; by++) {
      for(int bx=0; bx<blockScheme[0]; bx++) {
	int blk_idx = by * blockScheme[0] + bx;

	for(int i=0; i<nLow; i++) {	  
	  for(int k=0; k<blk_size; k++) {
	    block_data_temp[blk_idx][k] += block_coef[blk_idx][j*nLow + i] * block_data_ortho[blk_idx][blk_size * i + k];
	  }
	}	
      }
    }

    for(int by=0; by<blockScheme[1]; by++) {
      for(int bx=0; bx<blockScheme[0]; bx++) {
	int blk_idx = by * blockScheme[0] + bx;

	// Location of the start of the desired block (x runs fastest)
	int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));

	// Copy data to the jth eigenvector
	for(int ny=0; ny<y_block_size; ny++) {
	  for(int nx=0; nx<x_block_size; nx++) {
	    for(int mu=0; mu<2; mu++) {
	      
	      int loc_idx = 2*(nx + x_block_size * ny) + mu;      // Local index
	      int glo_idx = blk_offset + 2*((ny * Nx) + nx) + mu; // Global index
	      kSpace[j]->data[glo_idx] = block_data_temp[blk_idx][loc_idx];
	    }
	  }
	}
      }
    }
  }
} 
