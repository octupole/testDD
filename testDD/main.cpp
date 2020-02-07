
//============================================================================
// Name        : testMPI.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
#undef HAVE_MPI
#if HAVE_MPI
#include "../Parallel/MPI.h"
#include <MPIconfig.hpp>
#else
namespace Parallel{
static uint CARTDIRS{0};
struct MPI{
	void CartInit(){};
	MPI & gWorld(){return *this;}
	uint rank(){return 0;}
	uint size(){return 1;}
};
}
#endif
#include <map>
using namespace std;
template<class T, uint N>
struct helper
{
	static constexpr T pow(const T x){
		return helper<T, N-1>::pow(x) * x;
	}
};

template<class T>
struct helper<T, 0>
{
	static constexpr T pow(const T x){
		return 1;
	}
};
template<uint N, class T>
T constexpr pow(T const x)
{
	return helper<T, N>::pow(x);
}
constexpr double swrs(double & rout,double & rtol,double x){
	return 1-pow<2>((x-rout)/rtol)*(3-2*(x-rout)/rtol);
}
constexpr double dswrs(double & rout, double & rtol,double x){
	return -(2*(x-rout)/rtol)*(3-2*(x-rout)/rtol)/rtol+2.0*pow<2>((x-rout)/rtol)/rtol;
}

struct smoothstep{
	static double rout;
	static double rtol;

	static constexpr double sw(double x){
		return swrs(rout,rtol,x);
	}
	static constexpr double dsw(double x){
		return dswrs(rout,rtol,x);
	}
};
double smoothstep::rout=9;
double smoothstep::rtol=1;

constexpr double qp{0.3275911},a1{0.2548296},a2{-0.28449674},a3{1.4214137},a4{-1.453152},a5{1.0614054},twrtpi{2/sqrt(M_PI)};
constexpr double alpha_e(double & alpha){
	return alpha;
};
constexpr double alphar_e(double & alpha,double x){
	return alpha*x;
};

constexpr double qt_e(double & alpha,double x){
	return 1/(1+qp*alpha*x);
};
constexpr double expcst_e(double & alpha, double x){
	return exp(-pow<2>(alpha*x));
}

struct almost_erfc{
	static double alpha;
	static constexpr double alphar(double x){
		return alphar_e(alpha,x);
	};
	static constexpr double qt(double x){
		return qt_e(alpha,x);
	};
	static constexpr double expcst(double x){
		return expcst_e(alpha,x);
	}
	static constexpr double erfcst(double x){
		return ((((a5*qt(x)+a4)*qt(x)+a3)*qt(x)+a2)*qt(x)+a1)*qt(x)*expcst(x);
	}
	static constexpr double erfc_ri(double x){
		return ((((a5*qt(x)+a4)*qt(x)+a3)*qt(x)+a2)*qt(x)+a1)*qt(x)*expcst(x)/x;
	}

	static constexpr double derfcst(double x){
		return -twrtpi*alphar(1)*expcst(x);
	}
	static constexpr double derfc_ri(double x){
		return (x*derfcst(x)-erfcst(x))/pow<2>(x);
	}

};
double almost_erfc::alpha=0.3;

int main() {

//	Parallel::MPI my(10.0,100,100,100);
	Parallel::MPI my;
//	my.PrintInfo();
	my.CartInit();
	map<int,int> outbuf;
	for(int o{0};o<Parallel::CARTDIRS;o++){
		auto Cart=my.gWorld();
		if(o%2 == 0) outbuf[o]=Cart.rank();
		else outbuf[o]=-Cart.rank();
	}
	double rcut{10.0};
	smoothstep::rout=8.5;
	smoothstep::rtol=rcut-smoothstep::rout;
	almost_erfc::alpha=0.35;
	double dx=0.05;
	int stop=10/dx;
	for(size_t o{1};o<=stop;o++){
		double x=o*dx;
		auto der=-almost_erfc::erfcst(x)/(x*x)+almost_erfc::derfcst(x)/x;
		cout << x<< " " << almost_erfc::erfc_ri(x) << " " << almost_erfc::derfc_ri(x) <<endl;;
	}
//	my.CartSend(outbuf);

//	my.PrintInfo();
	return 0;
}
