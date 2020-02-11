/*
 * fftw3mpi.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: marchi
 */

#include "fftw3mpi.h"

namespace Parallel {


fftw3mpi::fftw3mpi(int nx, int ny, int nz) {
	int isInitialized;
	MPI_Initialized(&isInitialized);
	try{
	if(!isInitialized) throw string("MPI must be intialized before calling pfft!");
	}catch(const string & s){
		cout << s << endl;
		MPI_Finalize();
	}
	pfft_init();
	NN[XX]=nx;NN[YY]=ny;NN[ZZ]=nz;
}
fftw3mpi::fftw3mpi() {
	int isInitialized;
	MPI_Initialized(&isInitialized);
	try{
	if(!isInitialized) throw string("MPI must be intialized before calling pfft!");
	}catch(const string & s){
		cout << s << endl;
		MPI_Finalize();
	}
	pfft_init();
}
bool fftw3mpi::isInitialized(){
	bool ok=true;
	for(size_t o{0};o<DIM;o++){
		ok&=np[o]>0;
		ok&=NN[o]>0;
	}
	return ok;
}
void fftw3mpi::setDimFFT(int nx,int ny, int nz){
	NN[XX]=nx;NN[YY]=ny;NN[ZZ]=nz;
}
void fftw3mpi::getMemory(){
	try{
		if(!isInitialized()) throw string("fftw3mpi not initialized");
	}catch(const string & s){cout << s <<endl;MPI_Finalize();exit(0);}

	auto alloc_local=pfft_local_size_dft_r2c_3d(NN, *comm_cart_3dx, PFFT_TRANSPOSED_NONE,
	      local_ni, local_i_start, local_no, local_o_start);
	in.Allocate(2*alloc_local);
	out.Allocate(alloc_local);
}
void fftw3mpi::create3DMesh(int nx, int ny, int nz){
	np[XX]=nx;
	np[YY]=ny;
	np[ZZ]=nz;
	int npp[DIM]={nx,ny,nz};
	comm_cart_3dx=new MPI_Comm;

	try{
	 if(pfft_create_procmesh(DIM, MPI_COMM_WORLD,npp, comm_cart_3dx)) throw string("Error: This test file only works with ");
	}catch(const string & s){
		cout << s << np[XX]*np[YY]*np[ZZ] << " processes "<<endl;
		MPI_Finalize();exit(0);
	}
}
void fftw3mpi::create3DMesh(int n){
	vector<int> np0=findSize3D(n);
	np[XX]=np0[XX];
	np[YY]=np0[YY];
	np[ZZ]=np0[ZZ];

	int npp[DIM]={np0[XX],np0[YY],np0[ZZ]};
	comm_cart_3dx=new MPI_Comm;
	try{
	 if(pfft_create_procmesh(DIM,  MPI_COMM_WORLD, npp, comm_cart_3dx)) throw string("Error: This test file only works with ");
	}catch(const string & s){
		cout << s << np[XX]*np[YY]*np[ZZ] << " processes "<<endl;
		MPI_Finalize();exit(0);
	}
}
fftw3mpi::rcfft3d_mpi::rcfft3d_mpi(fftw3mpi & x){
	try{
		if(!x.isInitialized()) throw string("fftw3mpi not initialized");
	}catch(const string & s){cout << s <<endl;MPI_Finalize();exit(0);}
	myfftx=&x;
	double * in=&x.in[0];
	fftw_complex * out=(fftw_complex *) &x.out[0];

	rcfft3d_mpi::plan_dft=pfft_plan_dft_r2c_3d(x.NN,in, out, x.comm_cart_3d(),
			PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
}
fftw3mpi::crfft3d_mpi::crfft3d_mpi(fftw3mpi & x){
	myfftx=&x;
	double * in=&x.in[0];
	fftw_complex * out=(fftw_complex *) &x.out[0];

	crfft3d_mpi::plan_dft=pfft_plan_dft_c2r_3d(x.NN, out, in, x.comm_cart_3d(),
			PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
}
void fftw3mpi::rcfft3d_mpi::fft(array3<double> & A, array3<Complex> & B){
	double * in=&A[0][0][0];
	pfft_complex * out=(pfft_complex *)&B[0][0][0];
	pfft_execute_dft_r2c(rcfft3d_mpi::plan_dft,in,out);
}

void fftw3mpi::crfft3d_mpi::fft(array3<Complex> & B,array3<double> & A){
	double * in=&A[0][0][0];
	pfft_complex * out=(pfft_complex *)&B[0][0][0];
	pfft_execute_dft_c2r(crfft3d_mpi::plan_dft,out,in);
}
void fftw3mpi::crfft3d_mpi::fftnormalize(array3<Complex> & B,array3<double> & A){
	double * in=&A[0][0][0];
	pfft_complex * out=(pfft_complex *)&B[0][0][0];
	pfft_execute_dft_c2r(crfft3d_mpi::plan_dft,out,in);
//	for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++){
//	    in[l] /= (NN[XX]*NN[YY]*NN[ZZ]);
//	}

}

fftw3mpi::~fftw3mpi() {}

} /* namespace Parallel */
