#include "utils.h"
#include "io.h"

void writeGauge(field<Complex> *gauge, string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge->read(x,y,mu)) << endl;
  
  outPutFile.close();
  return;  
}

void readGauge(field<Complex> *gauge, string name)
{
  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	getline(inPutFile, val);
	gauge->write(x, y, mu, polar(1.0, stod(val)));
      }
    }
  }
  
  return;
}

#include "hdf5.h"

void hdf5Example() {

  hid_t       file_id, dataset_id, dataspace_id; // identifiers 
  herr_t      status;
  int         i, j;//, int_data[4][6], read_int_data[4][6];
  hsize_t     dims[2];
  short       short_data[4][6], read_short_data[4][6];
  
  // Initialize some datasets
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 6; j++) {
      short_data[i][j] = i * 6 + j + 1;
    }
  }

  // Create a new file using default properties. 
  file_id = H5Fcreate("short.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  // Create the data space for the dataset. 
  dims[0] = 4;
  dims[1] = 6;
  dataspace_id = H5Screate_simple(2, dims, NULL);

  // Create the dataset. 
  dataset_id = H5Dcreate2(file_id, "/short_data", H5T_STD_B16LE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  // Write the dataset. 
  status = H5Dwrite(dataset_id, H5T_STD_B16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    short_data);

  // End access to the dataset and release resources used by it. 
  status = H5Dclose(dataset_id);

  //------------------------------------------------------

  // Open an existing dataset. 
  dataset_id = H5Dopen2(file_id, "/short_data", H5P_DEFAULT);
  
  status = H5Dread(dataset_id, H5T_STD_B16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_short_data);
  
  // check data set
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 6; j++) {
      if(short_data[i][j] == read_short_data[i][j]) {
	//cout << "HDF5 woo hoo" << endl;
      } else if(short_data[i][j] != read_short_data[i][j]) {
	//cout << "HDF5 boo boo" << endl;
      }
    }
  }
  
  // Close the dataset. 
  status = H5Dclose(dataset_id);
  
  // Close the file. 
  status = H5Fclose(file_id);
}
