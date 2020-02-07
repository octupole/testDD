
//============================================================================
// Name        : testMPI.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, *in;
  pfft_complex *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_3d;

  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2; np[2] = 2;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_3d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_real(2 * alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft_r2c_3d(
      n, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft_c2r_3d(
      n, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_real(3, n, local_ni, local_i_start,
      in);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);

  /* clear the old input */
  pfft_clear_input_real(3, n, local_ni, local_i_start,
      in);

  /* execute parallel backward FFT */
  pfft_execute(plan_back);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_real(3, n, local_ni, local_i_start, in, comm_cart_3d);
  pfft_printf(comm_cart_3d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]);
  pfft_printf(comm_cart_3d, "maxerror = %6.2e;\n", err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_3d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
//
//#include <iostream>
//#include <vector>
//#include <cmath>
//#undef HAVE_MPI
//#if HAVE_MPI
//#include "../Parallel/MPI.h"
//#include <MPIconfig.hpp>
//#else
//namespace Parallel{
//static uint CARTDIRS{0};
//struct MPI{
//	void CartInit(){};
//	MPI & gWorld(){return *this;}
//	uint rank(){return 0;}
//	uint size(){return 1;}
//};
//}
//#endif
//#include <map>
//using namespace std;
//template<class T, uint N>
//struct helper
//{
//	static constexpr T pow(const T x){
//		return helper<T, N-1>::pow(x) * x;
//	}
//};
//
//template<class T>
//struct helper<T, 0>
//{
//	static constexpr T pow(const T x){
//		return 1;
//	}
//};
//template<uint N, class T>
//T constexpr pow(T const x)
//{
//	return helper<T, N>::pow(x);
//}
//constexpr double swrs(double & rout,double & rtol,double x){
//	return 1-pow<2>((x-rout)/rtol)*(3-2*(x-rout)/rtol);
//}
//constexpr double dswrs(double & rout, double & rtol,double x){
//	return -(2*(x-rout)/rtol)*(3-2*(x-rout)/rtol)/rtol+2.0*pow<2>((x-rout)/rtol)/rtol;
//}
//
//struct smoothstep{
//	static double rout;
//	static double rtol;
//
//	static constexpr double sw(double x){
//		return swrs(rout,rtol,x);
//	}
//	static constexpr double dsw(double x){
//		return dswrs(rout,rtol,x);
//	}
//};
//double smoothstep::rout=9;
//double smoothstep::rtol=1;
//
//constexpr double qp{0.3275911},a1{0.2548296},a2{-0.28449674},a3{1.4214137},a4{-1.453152},a5{1.0614054},twrtpi{2/sqrt(M_PI)};
//constexpr double alpha_e(double & alpha){
//	return alpha;
//};
//constexpr double alphar_e(double & alpha,double x){
//	return alpha*x;
//};
//
//constexpr double qt_e(double & alpha,double x){
//	return 1/(1+qp*alpha*x);
//};
//constexpr double expcst_e(double & alpha, double x){
//	return exp(-pow<2>(alpha*x));
//}
//
//struct almost_erfc{
//	static double alpha;
//	static constexpr double alphar(double x){
//		return alphar_e(alpha,x);
//	};
//	static constexpr double qt(double x){
//		return qt_e(alpha,x);
//	};
//	static constexpr double expcst(double x){
//		return expcst_e(alpha,x);
//	}
//	static constexpr double erfcst(double x){
//		return ((((a5*qt(x)+a4)*qt(x)+a3)*qt(x)+a2)*qt(x)+a1)*qt(x)*expcst(x);
//	}
//	static constexpr double erfc_ri(double x){
//		return ((((a5*qt(x)+a4)*qt(x)+a3)*qt(x)+a2)*qt(x)+a1)*qt(x)*expcst(x)/x;
//	}
//
//	static constexpr double derfcst(double x){
//		return -twrtpi*alphar(1)*expcst(x);
//	}
//	static constexpr double derfc_ri(double x){
//		return (x*derfcst(x)-erfcst(x))/pow<2>(x);
//	}
//};
//double almost_erfc::alpha=0.3;
//
//int main() {
//
////	Parallel::MPI my(10.0,100,100,100);
//	Parallel::MPI my;
////	my.PrintInfo();
//	my.CartInit();
//	map<int,int> outbuf;
//	for(int o{0};o<Parallel::CARTDIRS;o++){
//		auto Cart=my.gWorld();
//		if(o%2 == 0) outbuf[o]=Cart.rank();
//		else outbuf[o]=-Cart.rank();
//	}
//	double rcut{10.0};
//	smoothstep::rout=8.5;
//	smoothstep::rtol=rcut-smoothstep::rout;
//	almost_erfc::alpha=0.35;
//	double dx=0.05;
//	int stop=10/dx;
//	for(size_t o{1};o<=stop;o++){
//		double x=o*dx;
//		auto der=-almost_erfc::erfcst(x)/(x*x)+almost_erfc::derfcst(x)/x;
//		cout << x<< " " << almost_erfc::erfc_ri(x) << " " << almost_erfc::derfc_ri(x) <<endl;;
//	}
////	my.CartSend(outbuf);
//
////	my.PrintInfo();
//	return 0;
//}
