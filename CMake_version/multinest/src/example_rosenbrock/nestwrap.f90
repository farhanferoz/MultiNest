module nestwrapper

! Nested sampling includes

use Nested
use params
use like
   
implicit none
   
contains

!-----*-----------------------------------------------------------------

subroutine nest_Sample
	
	implicit none
	
   	integer nclusters				! total number of clusters found
	integer context
   	integer maxNode 				! variables used by the posterior routine
   
   
   	! calling MultiNest
	
   	call nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,nest_nPar, &
   	nest_nClsPar,nest_maxModes,nest_updInt,nest_Ztol,nest_root,nest_rseed,nest_pWrap, &
   	nest_fb,nest_resume,nest_outfile,nest_initMPI,nest_logZero,nest_maxIter,getLogLike,dumper,context)

end subroutine nest_Sample

!-----*-----------------------------------------------------------------

! Wrapper around Likelihood Function

subroutine getLogLike(Cube,n_dim,nPar,lnew,context)
	
	implicit none
	
	! Input arguments
	integer n_dim 					! dimensionality (total number of free parameters) of the problem
	integer nPar 					! total number of free plus derived parameters

	!Input/Output arguments
	double precision Cube(nPar) 			! on entry has the ndim parameters in unit-hypercube
	 						! on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
	 
	! Output arguments
	double precision lnew 				! loglikelihood
	integer context					! not needed, any additional information user wants to pass
	
   
   	
	!call your loglike function here 
	
   	call slikelihood(Cube,lnew)

end subroutine getLogLike

!-----*-----------------------------------------------------------------

! dumper, called after every updInt*10 iterations

subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSlogZ, logZerr, context)

	implicit none

	integer nSamples				! number of samples in posterior array
	integer nlive					! number of live points
	integer nPar					! number of parameters saved (physical plus derived)
	double precision, pointer :: physLive(:,:)	! array containing the last set of live points
	double precision, pointer :: posterior(:,:)	! array with the posterior distribution
	double precision, pointer :: paramConstr(:)	! array with mean, sigmas, maxlike & MAP parameters
	double precision maxLogLike			! max loglikelihood value
	double precision logZ				! log evidence value from the default (non-INS) mode
	double precision INSlogZ			! log evidence value from the INS mode
	double precision logZerr			! error on log evidence
	integer context					! not required by MultiNest, any additional information user wants to pass
	
end subroutine dumper

!-----*-----------------------------------------------------------------

end module nestwrapper
