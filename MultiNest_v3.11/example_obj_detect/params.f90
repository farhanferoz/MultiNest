! Include file for example MultiNest program obj_detect (see arXiv:0704.3704)

module params
implicit none

! Toy Model Parameters
	integer sdim
      	integer snpix !no. of pixels in the grid (generate a grid with snpix X snpix pixels)
      	integer snclstr !no. of objects to generate, should be a square no.
	parameter (snpix=200,snclstr=8)
      	logical autopos !position clusters uniformly in the field?
      	parameter (autopos=.false.)
	double precision sodata(snpix,snpix)
      	double precision spos(2*snclstr) !cluster position coordinates
      	data spos /0.7, 110.5, 68.2, 166.4, 75.3, 117.0, 78.6, &
      	12.6, 86.8, 41.6, 113.7, 43.1, 124.5, 54.2, 192.3, 150.2/
	double precision samp(snclstr) !amplitude of each object
      	data samp /0.71, 0.91, 0.62, 0.60, 0.63, 0.56, 0.60, 0.90/
	double precision ssig(snclstr) !width of each object
      	data ssig /5.34, 5.40, 5.66, 7.06, 8.02, 6.11, 9.61, 9.67/
	double precision snoise !gaussian noise rms
      	data snoise / 2. /
      	integer dseed!seed for creating mock data,-ve means take it from sys clock
	parameter(dseed=12)
      	double precision spriorran(4,2) !priors on the parameters
      	!uniform prior on x position
	data spriorran(1,1),spriorran(1,2) / 0. , 200. /
      	!uniform prior on y position
	data spriorran(2,1),spriorran(2,2) / 0. , 200. /
      	!uniform prior on amplitude
	data spriorran(3,1),spriorran(3,2) / 0. , 2. /
      	!uniform prior on sigma
	data spriorran(4,1),spriorran(4,2) / 3. , 12. /

! Parameters for MultiNest
	
      	!whether to do use Nested Importance Sampling
	logical nest_IS
 	parameter(nest_IS=.true.)
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.true.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.false.)
	
      	!max no. of live points
      	integer nest_nlive
	parameter(nest_nlive=1000)
      
      	!tot no. of parameters, should be sdim in most cases but if you need to
      	!store some additional parameters with the actual parameters then
      	!you need to pass them through the likelihood routine
	integer nest_nPar 
      
      	!seed for MultiNest, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
      	!evidence tolerance factor
      	double precision nest_tol 
      	parameter(nest_tol=0.5)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
      	parameter(nest_efr=0.8d0)
      
      	!root for saving posterior files
      	character*1000 nest_root
	parameter(nest_root='chains/2-')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=1000)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=10)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	parameter(nest_nClsPar=2)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	parameter(nest_resume=.true.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.true.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1d0))
      
      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	!parameters to wrap around (0 is F & non-zero T)
	integer nest_pWrap(4)
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================


end module params
