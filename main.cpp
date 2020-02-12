//#include <complex.h>
//#include <pfft.h>
//#include <iostream>
//
//using namespace std;
//int main(int argc, char **argv)
//{
//  int np[3];
//  ptrdiff_t n[3];
//  ptrdiff_t alloc_local,alloc_local0;
//  ptrdiff_t local_ni[3], local_i_start[3];
//  ptrdiff_t local_no[3], local_o_start[3];
//  double err, *in;
//  pfft_complex *out;
//  pfft_plan plan_forw=NULL, plan_back=NULL;
//  MPI_Comm comm_cart_3d;
//
//  /* Set size of FFT and process mesh */
//  n[0] = 4; n[1] = 4; n[2] = 4;
//  np[0] = 2; np[1] = 2; np[2] = 1;
//
//  /* Initialize MPI and PFFT */
//  MPI_Init(&argc, &argv);
//  pfft_init();
//
//  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
//  if( pfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
//    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]*np[2]);
//    MPI_Finalize();
//    return 1;
//  }
//
//  /* Get parameters of data distribution */
//  alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_3d, PFFT_TRANSPOSED_NONE,
//      local_ni, local_i_start, local_no, local_o_start);
//  int XX{0},YY{1},ZZ{2};
//
//  /* Allocate memory */
//  in  = pfft_alloc_real(2 * alloc_local);
//  out = pfft_alloc_complex(alloc_local);
//
//  /* Plan parallel forward FFT */
//  plan_forw = pfft_plan_dft_r2c_3d(
//      n, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
//
//  /* Plan parallel backward FFT */
//  plan_back = pfft_plan_dft_c2r_3d(
//      n, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
//
//  /* Initialize input with random numbers */
//  pfft_init_input_real(3, n, local_ni, local_i_start,
//      in);
//  int myrank, size;
//  ptrdiff_t m;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  ptrdiff_t *lis, *lni;
//
//  int nzp=n[2]/2+1;
//  lis = local_i_start; lni = local_ni;
//
//  /* Output results: here we want to see the data ordering of real and imaginary parts */
//  MPI_Barrier(MPI_COMM_WORLD);
//  for(int t=0; t<size; t++){
//	  if(myrank == t){
//		  printf("rank %d: R2C PFFT Input:\n", myrank);
//		  printf("rank %d: lis = [%td, %td, %td], lni = [%td, %td, %td]\n", myrank, lis[0], lis[1], lis[2], lni[0], lni[1], lni[2]);
//		  m=0;
//		  for(ptrdiff_t k0=lis[0]; k0<lis[0]+lni[0]; k0++)
//			  for(ptrdiff_t k1=lis[1]; k1<lis[1]+lni[1]; k1++)
//				  for(ptrdiff_t k2=lis[2]; k2<lis[2]+lni[2]; k2++, m++){
//					  int mm=n[1]*n[2]*k0+n[2]*k1+k2;
//					  printf("%d <--> %d  in[%td, %td, %td] = %.2f\n", m, mm, k0, k1, k2, in[m]);
//				  }
//		  fflush(stdout);
//	  }
//	  MPI_Barrier(MPI_COMM_WORLD);
//  }
//
//  /* execute parallel forward FFT */
//  pfft_execute(plan_forw);
//
//  /* clear the old input */
//  pfft_clear_input_real(3, n, local_ni, local_i_start,
//      in);
//
//
//    ptrdiff_t *los, *lno;
//    los = local_o_start; lno = local_no;
//
//    /* Output results: here we want to see the data ordering of real and imaginary parts */
//    MPI_Barrier(MPI_COMM_WORLD);
//    for(int t=0; t<size; t++){
//      if(myrank == t){
//        printf("rank %d: R2C PFFT Output:\n", myrank);
//        printf("rank %d: los = [%td, %td, %td], lno = [%td, %td, %td]\n", myrank, los[0], los[1], los[2], lno[0], lno[1], lno[2]);
//        m=0;
//        for(ptrdiff_t k0=los[0]; k0<los[0]+lno[0]; k0++)
//          for(ptrdiff_t k1=los[1]; k1<los[1]+lno[1]; k1++)
//            for(ptrdiff_t k2=los[2]; k2<los[2]+lno[2]; k2++, m++){
//            	int mm=ny*nzp*k0+nzp*k1+k2;
//            	printf("%d <--> %d  out[%td, %td, %td] = %.2f + I * %.2f\n", m, mm, k0, k1, k2, out[m][0], out[m][1]);
//            }
//        fflush(stdout);
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
//  /* execute parallel backward FFT */
//  pfft_execute(plan_back);
//
//  /* Scale data */
//  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
//    in[l] /= (n[0]*n[1]*n[2]);
//
//  /* Print error of back transformed data */
//  MPI_Barrier(MPI_COMM_WORLD);
//  err = pfft_check_output_real(3, n, local_ni, local_i_start, in, comm_cart_3d);
//  pfft_printf(comm_cart_3d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]);
//  pfft_printf(comm_cart_3d, "maxerror = %6.2e;\n", err);
//
//  /* free mem and finalize */
//  pfft_destroy_plan(plan_forw);
//  pfft_destroy_plan(plan_back);
//  MPI_Comm_free(&comm_cart_3d);
//  pfft_free(in); pfft_free(out);
//  MPI_Finalize();
//  return 0;
//}

#include <complex.h>
#include <pfft.h>
#include <iostream>
#include <random>

#include "../Parallel/fftw3mpi.h"
#include "../Tools/Array.h"
#include "../Tools/Ftypedefs.h"
using std::cout;
using std::endl;
using namespace Array;
using Complex=std::complex<double>;
int main(int argc, char **argv)
{

  int np[3];
  ptrdiff_t n[3];
  /* Set size of FFT */
  n[XX]=128;n[YY]=128;n[ZZ]=128;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  int psize, myrank;
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  Parallel::fftw3mpi myFFT(n[XX],n[YY],n[ZZ]);
  myFFT.create3DMesh(psize);
  myFFT.getMemory();

  auto lis = myFFT.loc_i_start();
  auto lni =myFFT.loc_ni();
  auto los = myFFT.loc_o_start();
  auto lno = myFFT.loc_no();

  int nzp=n[ZZ]/2+1;
  array3<double> A(n[XX],n[YY],n[ZZ]);
  array3<double> C(n[XX],n[YY],n[ZZ]);
  array3<Complex> B(n[XX],n[YY],nzp);
  myFFT.init(A);

  C=A;

  Parallel::fftw3mpi::rcfft3d_mpi forward(myFFT);
  Parallel::fftw3mpi::crfft3d_mpi backward(myFFT);

  /* execute parallel forward FFT */
  double t1=MPI_Wtime();

  forward.fft(A,B);


  backward.fftnormalize(B,A);

  /* Print error of back transformed data */

  double t2=MPI_Wtime();
//  if(!myrank){
	  cout << "Timing "<< t2-t1  <<endl;
//  }
  double err=-1.0e12;
  size_t m=0;
  for(int k0=lis[0]; k0<lis[0]+lni[0]; k0++)
    for(int k1=lis[1]; k1<lis[1]+lni[1]; k1++)
      for(int k2=lis[2]; k2<lis[2]+lni[2]; k2++,m++){
    	  double dx=A[k0][k1][k2]-C[k0][k1][k2];
    	  if(fabs(dx)>err)err=fabs(dx);
      }
  fflush(stdout);
  vector<double> myerr(psize);
  MPI_Gather(&err,1,MPI_DOUBLE, &myerr[0],1,MPI_DOUBLE,0, MPI_COMM_WORLD);
  err=-1.0e12;
  if(myrank == 0) {
	  for(auto it: myerr)
		  if(it > err) err=it;
	  cout << "Error "<< err<<endl;
  }
  MPI_Finalize();
  return 0;
}
