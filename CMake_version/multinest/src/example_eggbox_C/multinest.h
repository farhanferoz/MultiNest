#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

#ifndef MULTINEST_H
#define MULTINEST_H

#ifdef __cplusplus

/***************************************** C++ Interface to MultiNest **************************************************/

#include <cstring>

namespace nested
{

	// map the Fortran 90 entry points of libnest3.a to C++ functions

	// module nested, function nestRun maps to nested::run

	// the pass-by-reference nature of most of the Fortran is translated away
	// *apart* from the callbacks. The provided call back functions must still accept 
	// references rather than values. There is also some confusion as to the type
	// of the first argument of LogLike. 
	// Should it be a double * or an farray<double, 1> *? The former seems to 
	// work and is simpler.

	// This structure is reverse engineered from looking 
	// at gfortran stack traces. It is probably wrong
	
	template<typename type, int ndims> class farray_traits;
	
	template<> class farray_traits<double, 1> { public: static const int id = 537; };
	template<> class farray_traits<double, 2> { public: static const int id = 538; };
	template<> class farray_traits<int, 1> { public: static const int id = 265; };
	template<> class farray_traits<int, 2> { public: static const int id = 266; };

	// the extra data for f90 that defines how arrays are arranged.
	template<typename T, int ndim> class farray
	{
		public:
			farray(T *_data, int w, int h = 0) : data(_data), offset(0), type(farray_traits<T, ndim>::id), 
			x_stride(1), x_lbound(1), x_ubound(w), y_stride(w), y_lbound(1), y_ubound(h) {};
			
			T *data;
			int offset;
			int type;
			int x_stride, x_lbound, x_ubound;
			int y_stride, y_lbound, y_ubound;
	};
	
	extern "C" {
		void NESTRUN(int &IS, int &mmodal, int &ceff, int &nlive, double &tol, double &efr, int &ndims,
			int &nPar, int &nClsPar, int &maxModes, int &updInt, double &Ztol, char *root, int &seed,
			int *pWrap, int &fb, int &resume, int &outfile, int &initMPI, double &logZero, int &maxiter,
			void (*Loglike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
			void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *),
			void *context, int &root_len);
	}

	static void run(bool IS, bool mmodal, bool ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, int maxModes,
		int updInt, double Ztol, const std::string & root, int seed, int *pWrap, bool fb, bool resume, bool outfile, 
		bool initMPI, double logZero, int maxiter, void (*LogLike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
		void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *), void *context)
	{
		char t_root[1000];
		std::fill(t_root, t_root + 1000, ' ');
		snprintf(t_root, 999, "%s", root.c_str());
		int root_len = strlen(t_root);
		t_root[strlen(t_root)] = ' ';
	
		int t_fb = fb;
		int t_resume = resume;
		int t_outfile = outfile;
		int t_initMPI = initMPI;
		int t_mmodal = mmodal;
		int t_IS = IS;
		int t_ceff = ceff;
		
		NESTRUN(t_IS, t_mmodal, t_ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, t_root, seed, pWrap, t_fb, 
		t_resume, t_outfile, t_initMPI, logZero, maxiter, LogLike, dumper, context, root_len);
	}	
}

/***********************************************************************************************************************/

#else // ifdef __cplusplus

/***************************************** C Interface to MultiNest **************************************************/

extern void NESTRUN(int *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, 
char *, int *, int *, int *, int *, int *, int *, double *, int *, void (*Loglike)(double *, int *, int *, 
double *, void *), void (*dumper)(int *, int *, int *, double **, double **, double **, double *, 
double *, double *, double *, void *), void *context);

void run(int IS, int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, 
int maxModes, int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, 
int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, void *), 
void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, double *, void *), 
void *context)
{
	int i;
	for (i = strlen(root); i < 1000; i++) root[i] = ' ';

        NESTRUN(&IS, &mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/

#endif // ifdef __cplusplus

#endif // MULTINEST_H
