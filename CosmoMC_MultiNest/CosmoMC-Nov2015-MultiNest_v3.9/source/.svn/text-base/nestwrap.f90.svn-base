! The wrapper for MultiNest 3.5, CosmoMC Mar 2014
! Author: Farhan Feroz & Will Handley (wh260@cam.ac.uk) Mar 2014

module nestwrap

    use settings
    use BaseParameters
    use GeneralSetup
    use nested, DUMMY1 => MPI_COMM_WORLD
    use MpiUtils

    implicit none

    Type :: Prior
        ! This derived type has the aim of simplifying the dealings with the
        ! parameter indices.
        ! Normal cosmomc has an integer array params_used(:)
        ! These integers refer to the array indices of parameters to be varied.
        ! These parameters are stored in P within the derived type ParamSet.
        ! Thus CurParams%P(params_used(i)) refers to the ith used parameter
        ! stored in CurParams%P(:)

        ! This type aims to split up params_used(:) into the separate types of
        ! priors
        ! Imagine there are three different type of priors, A,B,C (for example
        ! Uniform, Log-Uniform and Gaussian) on 20 variables. 
        !
        ! We split the 20 Cube parameters into the three distict sets, the
        ! indices of the cube variable stored in their cube_pars(:) array, and end
        ! position of the previous region stored in the variable cube_pos
        !
        ! Further, the corresponding indices of the cosmological parameters are
        ! stored in the cube_pos(:) array.
        !
        ! A%cube_pos   B%cube_pos        C%cube_pos
        ! |                |                 |
        ! v                v                 v
        ! | 1  2  3  4  5  6  7  8  9 10 11 12 13 14 16 17 18 19 20|   Cube(:)
        ! |<--A%cube_pars-->|<--B%cube_pars-->|<---C%cube_pars---->|
        !
        !
        ! | 1 2 3 4 16 19 22 23 24 25 26 27 28 29 30 31 32 33 34 35| params_used(:)
        ! |<-A%cosmo_pars-->|<-B%cosmo_pars-->|<---C%cosmo_pars--->|
        !
        !
        ! Thus, if you want to transfer the ith parameter of type A from the
        ! cube to cosmological parameters via a function f, one can do that
        ! with:
        ! CurParams%P(A%cosmo_pars(i))   =   f(   Cube%P(A%cube_pars(i))  )

        integer :: num_var   = 0 ! The number of variables stored
        integer :: cube_pos  = 0 !
        integer, dimension(:), allocatable   :: cosmo_pars
        integer, dimension(:), allocatable   :: cube_pars
        contains
        procedure :: initialize_prior
        procedure :: add_point
        procedure :: set_cube_pars
    end Type


    !nested sampling parameters
    logical nest_IS !do importance nested sampling?
    logical nest_mmodal !multiple modes expected?
    logical nest_ceff !run in constant efficiency mode?
    integer nest_nlive !no. of live points
    double precision nest_tol !evidence tolerance factor
    double precision nest_ef !sampling efficiency
    integer nest_ndims !dimensionality
    integer nest_totPar !tot no. of parameters to be saved along with the sampling parameters
    integer nest_nCdims !total no. of parameters on which clustering should be done
    integer maxClst !max modes expected, for memory allocation
    integer nest_updInt !no. of iterations after which to update on the progress
    double precision nest_Ztol !lowest local evidence for which samples to produce
    character*100 nest_root !root for saving posterior files
    integer seed !seed for nested sampler, -ve means take it from sys clock
    integer, dimension(:), allocatable :: nest_pWrap ! Cyclic parameters
    logical nest_fb !feedback on the sampling progress?
    logical nest_resume !resume from a previous run?
    logical nest_outfile !Produce MultiNest output?
    logical initMPI !Need to initialise MPI?
    double precision, parameter :: nest_logZero = -huge(1.d0)*epsilon(1.d0) !likelihood zero
    double precision, parameter :: nest_logUnphysical = -1.d30 ! The likelihood value equal and beneath which  
                                                               ! points are discarded as being 'unphysical'
    integer nest_maxIter !Maximum number of iterations of the algorithm

    integer nest_num_derived ! Number of derived parameters
    real(mcp), allocatable :: derived(:) ! Pointer used to calculate derived parameters

    !retain CurParams throughout, for the code to be run in the only fast params mode
    Type(ParamSet) CurParams

    Type(Prior) Uniform    ! index information for Uniform priors 
    Type(Prior) Gaussian   ! index information for Gaussian priors 
    Type(Prior) MvG        ! index information for MultiVariate Gaussian priors
    Type(Prior) Ellipsoid  ! index information for Uniform Ellipsoidal priors
    Type(Prior) Sorted     ! index information for Sorted Uniform priors


    ! File names for Multivariate gaussian files
    character(LEN=:), allocatable :: file_MvG_mean
    character(LEN=:), allocatable :: file_MvG_cov
    ! Store of MvG parameter names from file
    character(LEN=:), allocatable :: MvG_parameters
    ! Matrix for holding covariances and means
    double precision, dimension(:,:), allocatable   :: MvG_prior_matrix
    double precision MvG_stretch ! The stretch factor for the covariance matrix

    ! File names for Ellipsoidal uniform files
    character(LEN=:), allocatable :: file_Ellipsoid_mean
    character(LEN=:), allocatable :: file_Ellipsoid_cov
    ! Store of Ellipsoid parameter names from file
    character(LEN=:), allocatable :: Ellipsoid_parameters
    ! Matrix for holding covariances and means
    double precision, dimension(:,:), allocatable   :: Ellipsoid_prior_matrix
    double precision Ellipsoid_stretch ! The stretch factor for the covariance matrix 


    character(LEN=:), allocatable :: sorted_uniform_parameters ! Store of Sorted parameter names from file 
    double precision sorted_uniform_max ! The maximum for the sorted uniform parameters
    double precision sorted_uniform_min ! The minimum for the sorted uniform parameters



    real(mcp), allocatable ::  nest_GaussPriors_std(:)  ! Storage for the cosmoMC gaussian priors
    real(mcp), allocatable ::  nest_GaussPriors_mean(:) ! Storage for the cosmoMC gaussian priors

    integer prior_boundary ! keeps track of boundary of where to get parameters from Cube(:)
    double precision prior_vol
    character(LEN=:), allocatable :: prior_vol_file

    integer size_derived

    logical cyclic_boundaries ! whether or not to set all boundaries as cyclic


    ! Code parameters:
    !
    ! This number determines how we define what the 'volume' of the prior space
    ! a gaussian curve occupies. Technically a gaussian goes over the entire
    ! space, but if we only demand that we measure it by a percentage of it's
    ! probability massy, we find we can set the parameters as:
    !   std_fac = 1.d0    =>    66%
    !   std_fac = 2.d0    =>    95%
    !   std_fac = 3.d0    =>    99%
    ! However, the 'best' prior width is in fact sqrt(2*pi). This is due to the
    ! fact that it best quantifies how the prior volume enters the evidence
    ! calculation. For more details see pp 348-350 of Information Theory,
    ! Inference and Learning Algorithms (MacKay)
    ! http://wh260.user.srcf.net/Books/Information_Theory_Inference_and_Learning_Algorithms.pdf
    double precision, parameter :: std_fac = 2.50662827463100050241576528481

    ! The fractional error in our prior volume calculation
    double precision, parameter :: prior_vol_frac_err = 1d-2


    contains

    !-------------------------------------------------------------------------

    subroutine read_multinest_inputs(Ini)

        implicit none
        Type(TSettingIni) :: Ini

        !set up MultiNest parameters
        nest_IS = Ini%Read_Logical('do_INS',.true.)
        nest_mmodal = Ini%Read_Logical('multimodal',.false.)
        nest_ceff = Ini%Read_Logical('ceff',.false.)
        nest_nlive = Ini%Read_Int('nlive',400)
        nest_tol = Ini%Read_Double('tol',0.5d0)
        nest_ef = Ini%Read_Double('eff',0.3d0)
        maxClst = Ini%Read_Int('maxmodes',1)
        nest_updInt = Ini%Read_Int('updInt',50)
        nest_root=rootname

        sorted_uniform_parameters = Ini%Read_String_Default('sorted_uniform_parameters','')
        sorted_uniform_max = Ini%Read_Double('sorted_uniform_max',1.d0) 
        sorted_uniform_min = Ini%Read_Double('sorted_uniform_min',0.d0) 

        MvG_parameters      = Ini%Read_String_Default('multivariate_gaussian_parameters','')
        file_MvG_cov  = Ini%Read_String_Default('multivariate_gaussian_cov_mat'   ,'')
        file_MvG_mean = Ini%Read_String_Default('multivariate_gaussian_mean_vec'  ,'')
        MvG_stretch= Ini%Read_Double('multivariate_gaussian_stretch',1.0d0)

        Ellipsoid_parameters      = Ini%Read_String_Default('ellipsoidal_uniform_parameters','')
        file_Ellipsoid_cov  = Ini%Read_String_Default('ellipsoidal_uniform_cov_mat'   ,'')
        file_Ellipsoid_mean = Ini%Read_String_Default('ellipsoidal_uniform_mean_vec'  ,'')
        Ellipsoid_stretch= Ini%Read_Double('ellipsoidal_uniform_stretch',1.0d0)


        cyclic_boundaries = Ini%Read_Logical('cyclic_boundary',.false.)


        !set up the seed for MultiNest
        seed=-1 !take it from system clock

        !feedback for MultiNest
        if(FeedBack>0) then
            nest_fb=.true.
        else
            nest_fb=.false.
        end if

        ! Set the zero of the loglikelihood.
        nest_Ztol=-1.d90

        ! Produce standard MultiNest outputs (compatible with getdist)
        nest_outfile=.true.

        ! CosmoMC already initialises MPI in driver.f90, so MultiNest doesn't need to
        initMPI=.false.

        ! No maximum to the number of iterations of the algorithm
        nest_maxIter=-1


    end subroutine read_multinest_inputs
    !-------------------------------------------------------------------------

    subroutine setup_multinest

        implicit none

        double precision loglikelihood
        integer i
#ifdef MPI
        integer ierror
#endif

        !retain CurParams throughout, for the code to be run in the only fast params mode
        CurParams%P(:num_params) = BaseParams%center(:num_params)
        loglikelihood = -Setup%LikeCalculator%GetLogLike(CurParams)


        ! set dimensionality
        nest_ndims       = num_params_used

        ! find out the number of derived parameters
        call Setup%Config%Parameterization%CalcDerivedParams(CurParams%P,CurParams%Theory, derived) 
        call DataLikelihoods%addLikelihoodDerivedParams(CurParams%P, CurParams%Theory, derived, CurParams%Likelihoods, -loglikelihood)

        if (allocated(derived)) then
            size_derived=size(derived)
            deallocate(derived)
        else
            size_derived=0
        end if

        nest_num_derived = size_derived + DataLikelihoods%Count + DataLikelihoods%LikelihoodTypeIndices%Count

        !total number of parameters for multinest to carry around
        nest_totPar      = nest_ndims + nest_num_derived

        ! Give feedback on the above variables if desired
        if (Feedback > 1) then
            write(*,*) 'Dimensionality of param space = ', nest_ndims
            write(*,*) 'Number of derived parameters  = ', nest_num_derived
            write(*,*) 'Location of ln(like) in live points file, column = ', nest_totPar+1
        endif

        if (nest_mmodal) then
            ! If it's multimodal, we should perform modal separation on all of
            ! the parameters.
            nest_nCdims=nest_ndims
        else
            nest_nCdims=0
        endif

        ! Set up the array containing cyclic variable information
        allocate(nest_pWrap(nest_ndims))

        ! At Farhan's advice, a neat trick to speed up the initial convergence
        ! is to set all of the variables to be 'wraparound'.
        ! This means that if any of the ellipsoids are sampling outside of the
        ! prior volume, they are wrapped back in, and at least provide some
        ! useful information (with INS)
        if (cyclic_boundaries) then
            nest_pWrap=1
        else
            nest_pWrap=0
        endif




        ! Set the end point of the prior region of the cube to be 0
        prior_boundary = 0

        ! Zero the total volume of the priors
        prior_vol = 1



        !----------- Configure the Non-Separable Priors -----------------
        if (trim(MvG_parameters)/='') call setup_MvG_priors(prior_boundary,prior_vol)
        if (trim(Ellipsoid_parameters)/='') call setup_Ellipsoid_priors(prior_boundary,prior_vol)
        if (trim(sorted_uniform_parameters)/='') call setup_Sorted_priors(prior_boundary,prior_vol)

        !----------- Configure the Separable Priors ---------------------
        ! Remove cosmomc's native gaussian prior code
        call setup_Gaussian_priors(prior_boundary,prior_vol)

        ! set up the separable uniform and gaussian priors
        call setup_Separable_priors(prior_boundary,prior_vol)

        !----------- Calculate the modified volume ---------------------
        ! Since we're running in parallel anyway and this code takes a few
        ! seconds, we might as well parellise it (since it trivially does so
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
        call Trim_Volume(nest_ndims,nest_totPar,prior_vol)  
        ! Save the prior volume to file
        prior_vol_file=rootname//'.prior_vol'
        open(unit=5555,file=prior_vol_file)
        write(5555,*) prior_vol
        close(5555)

        ! Output information to command line
        if ( IsMainMPI() ) call write_priors_information(prior_vol)




    end subroutine setup_multinest

    !-------------------------------------------------------------------------


    subroutine write_priors_information(prior_vol)
        implicit none
        double precision prior_vol
        integer i                   !Counter for loop

        write(*,*) 'Priors information:'

        if (MvG%num_var /= 0) then
            write(*,'(" multivariate gaussian prior on ",1I2," parameters ")') MvG%num_var
            do i=1,MvG%num_var
                write(*,*) ' ',trim(BaseParams%NameMapping%NameAtIndex(MvG%cosmo_pars(i)))
            enddo
        endif

        if (Ellipsoid%num_var /= 0) then
            write(*,'(" ellipsoidal uniform prior on ",1I2," parameters ")') Ellipsoid%num_var
            do i=1,Ellipsoid%num_var
                write(*,*) ' ',trim(BaseParams%NameMapping%NameAtIndex(Ellipsoid%cosmo_pars(i)))
            enddo
        endif

        if (Sorted%num_var /= 0) then
            write(*,'(" sorted uniform priors on ",1I2," parameters ")') Sorted%num_var
            do i=1,Sorted%num_var
                write(*,*) ' ',trim(BaseParams%NameMapping%NameAtIndex(Sorted%cosmo_pars(i)))
            enddo
        endif

        if (Gaussian%num_var /= 0) then
            write(*,'(" separable gaussian priors on ",1I2," parameters ")') Gaussian%num_var
            do i=1,Gaussian%num_var
                write(*,*) ' ',trim(BaseParams%NameMapping%NameAtIndex(Gaussian%cosmo_pars(i)))
            enddo
        endif

        if (Uniform%num_var /= 0) then
            write(*,'(" separable uniform priors on ",1I2," parameters ")') Uniform%num_var
            do i=1,Uniform%num_var
                write(*,'(" ",A15,"  [",1f10.5,","1f10.5,"]")') &
                    trim(BaseParams%NameMapping%NameAtIndex(Uniform%cosmo_pars(i))),&
                    BaseParams%Pmin(Uniform%cosmo_pars(i)),&
                    BaseParams%Pmax(Uniform%cosmo_pars(i))
            enddo
        endif

        !write (*,'(a,1f7.4,a,1f5.2)') 'log(Prior Volume) (99% of gaussian):', log(prior_vol), '+/-', prior_vol_frac_err
        write (*,'(a,1f8.2,a,1f5.2)') 'log(Prior Volume) (99% of gaussian):', log(prior_vol), ' +/-', prior_vol_frac_err
    end subroutine write_priors_information



    subroutine getLogLikeNS(Cube,n_dim,nPar,loglikelihood,context)
        !
        ! Multinest is a uniform sampler on the unit hypercube, calculating both
        ! posteriors and the evidence for a likelihood defined on the hypercube.
        ! This function is Multinest's likelihood function.
        !
        ! Inputs:
        !   Cube(:) position in the multinest hypercube
        ! Outputs:
        !   loglikelihood    the new loglikelihood
        !   Cube(:) physical + derived parameters
        !
        ! This function's layout is thus:
        !
        ! 1) Transform each of the Multinest parameters [Cube(j)] from the range
        !    [0,1] to physical parameters. This tranformation is dependent on
        !    the priors (usually uniform, but don't have to be). The code for
        !    transforming parameters for various priors is in
        !    multinest/priors.f90
        !
        ! 2) Calculate the log(Likelihood) using the standard cosmomc code
        !    acting on these transformed parameters.
        !
        ! 3) Calculate the derived parameters from the physical parameters
        !
        ! 4) Hand the physical + derived parameters back to multinest (via
        !    Cube) for printing at the appropriate time. 


        implicit none
        ! Inputs:
        double precision Cube(nPar) !Current position in the multinest hypercube
        integer n_dim               !Dimensionality of the multinest space
        integer nPar                !Total number of parameters (physical + derived)
        ! Output:
        double precision loglikelihood       !New likelihood to be calculated
        ! Other:
        integer context             !Any additional information one would like to pass through MultiNest
        integer i                   !Counter for loop

        ! Local variables
        Type(ParamSet) TrialParams  !Temp storage for Params


        ! Set TrialParams to have the default CurParams
        TrialParams=CurParams
        ! setting cosparams to the initial value (as default)
        TrialParams%P(:num_params) = BaseParams%center(:num_params)

        ! 1) Transform the Cube parameters into cosmological parameters
        call Transform_Cube_to_Cosmo(Cube,TrialParams%P,nPar)

        ! 2) Calculate the Likelihood
        ! ---------------------------
        !
        ! Calculate the likelihood 

        if (any( TrialParams%P(params_used(:num_params_used)) < BaseParams%Pmin(params_used(:num_params_used))) .or. any( TrialParams%P(params_used(:num_params_used)) > BaseParams%Pmax(params_used(:num_params_used))) ) then
            ! Some choices of priors (MvG,Ellipsoid) may call the parameter
            ! outside of their physically allowed region. Since cosmoMC comes
            ! with default prior widths in Pmax and Pmin (contained in batch1/)
            ! we can use these as cutoffs for any other distribution.

            ! Set the loglike to the zero likelihood, causing MultiNest to
            ! discard the point and draw a new one
            loglikelihood=nest_logZero
        else

            ! Calculate the likelihood
            loglikelihood = -Setup%LikeCalculator%GetLogLike(TrialParams)


            ! 3) Calculate the derived parameters
            ! -----------------------------------
            ! This code has been highjacked from the function:
            !   TCalculationAtParamPoint_WriteParams
            ! found in source/GeneralTypes.f90
            !
            ! It serves to calculate the derived parameters and the individual
            ! likelihoods
            !
            call Setup%Config%Parameterization%CalcDerivedParams(TrialParams%P,TrialParams%Theory, derived) 
            call DataLikelihoods%addLikelihoodDerivedParams(TrialParams%P, TrialParams%Theory, derived)

            ! 4) Hand the physical and derived parameters back to multinest
            ! -------------------------------------------------------------
            ! Now we pass all of the parameters (actual + derived) to the Cube,
            ! which then goes through the rest of the MultiNest file to eventually
            ! be printed.

            ! Physical parameters (being varied)
            Cube(:num_params_used) = TrialParams%P(params_used(:num_params_used))

            ! Derived parameters
            if (size_derived>0) Cube(num_params_used+1:num_params_used+size_derived) =  derived

            ! Likelihoods
            Cube(num_params_used+size_derived+1:num_params_used+size_derived+DataLikelihoods%Count) = &
                  TrialParams%Likelihoods(1:DataLikelihoods%Count)*2

            do i=1, DataLikelihoods%LikelihoodTypeIndices%Count
                Cube(num_params_used+size_derived+DataLikelihoods%Count+i) = &
                    sum(TrialParams%Likelihoods(DataLikelihoods%LikelihoodTypeIndices%Item(i)))*2
            end do

            ! 5) Clean up
            ! ---------------------

            ! free up the memory which was allocated in the CalcDerivedParams
            ! function. (Need to do this since we call getloglikeNS a lot)
            !
            deallocate(derived)
        endif


        ! Set it to our zero if it's too low
        if(loglikelihood<=nest_logUnphysical) loglikelihood=nest_logZero








        ! Calling this subroutine frees the memory of the variable 'TrialParams'
        ! choosing false ensures that TrialParams is freed, and not CurParams
        !
        call CurParams%AcceptReject(TrialParams,.false.)



    end subroutine getLogLikeNS





    !-------------------------------------------------------------------------
    subroutine Transform_Cube_to_Cosmo(Cube,P,nPar)
        ! Transform Multinest Hypercube to Cosmological Parameters
        ! -----------------------------------------------------------
        ! The cube parameters are first split into the types of prior 
        !  A) Uniform priors                    (Uniform)    
        !  B) Gaussian Priors                   (Gaussian)   
        !  C) Multivariate Gaussian Priors      (MvG)        
        !  D) Ellipsoidal Uniform Priors        (MvG)        
        !  E) Sorted Uniform Priors             (Sorted)        
        !
        ! The information about where parameters are stored is kept in the
        ! 'Prior' types indicated in parenthesis above
        implicit none

        double precision, dimension(nPar), intent(in)  :: Cube(nPar)
        integer                          , intent(in)  :: nPar
        double precision, dimension(nPar), intent(out) :: P(nPar)

        double precision temp_vec(nPar)
        integer i


        ! A) Uniform priors
        !------------------
        !
        do i=1,Uniform%num_var

            ! Re scale the variable using the Uniform Prior re-scaling
            P(Uniform%cosmo_pars(i))            &
                =                                             &
                UniformPrior( Cube(Uniform%cube_pars(i) ),  &
                BaseParams%Pmin(Uniform%cosmo_pars(i)),         &
                BaseParams%Pmax(Uniform%cosmo_pars(i))          &
                )                                             
        enddo


        ! B) Standard Gaussian Priors
        !----------------------------
        !
        do i=1,Gaussian%num_var

            ! Re scale the variable using the Gaussian Prior re-scaling
            P(Gaussian%cosmo_pars(i))                   &
                =                                                     &
                GaussianPrior(  Cube(Gaussian%cube_pars(i)),        &
                nest_GaussPriors_mean(Gaussian%cosmo_pars(i)),      &
                nest_GaussPriors_std (Gaussian%cosmo_pars(i))       &
                )
        enddo


        ! C) Multivariate Gaussian Priors
        !--------------------------------
        ! Transform the vector of gaussian prior distributed parameters from the
        ! MultiNest Hypercube to the physical space
        !       - defined by a covariance matrix across several parameters
        !       - this will most likely have come from a previous cosmological
        !         run (e.g. WMAP)
        ! The variable MvG_prior_matrix has the means in its last column, and
        ! the covariances in the rest of the matrix
        !
        if (MvG%num_var /=0 ) then
            temp_vec(1:MvG%num_var) = Cube(MvG%cube_pars)
            call MultivariateGaussianPrior( temp_vec(1:MvG%num_var), &
                MvG_prior_matrix(MvG%num_var+1,:) ,                  &
                MvG_prior_matrix(1:MvG%num_var,1:MvG%num_var)) 

            P(MvG%cosmo_pars) = temp_vec(1:MvG%num_var) 
        endif

        ! D) Ellipsoidal Uniform Priors
        !--------------------------------
        ! Transform the vector of Ellipsoidal uniform prior distributed parameters from the
        ! MultiNest Hypercube to the physical space
        !       - defined by a covariance matrix across several parameters
        !       - this will most likely have come from a previous cosmological
        !         run (e.g. WMAP)
        ! The variable Ellipsoid_prior_matrix has the means in its last column, and
        ! the covariances in the rest of the matrix
        !
        if (Ellipsoid%num_var /= 0 ) then
            temp_vec(1:Ellipsoid%num_var) = Cube(Ellipsoid%cube_pars)
            call EllipsoidalUniformPrior( temp_vec(1:Ellipsoid%num_var), &
                Ellipsoid_prior_matrix(Ellipsoid%num_var+1,:) ,                  &
                Ellipsoid_prior_matrix(1:Ellipsoid%num_var,1:Ellipsoid%num_var)) 

            P(Ellipsoid%cosmo_pars) = temp_vec(1:Ellipsoid%num_var) 
        endif

        ! E) Sorted Uniform Priors
        !--------------------------------
        ! Transform the vector of variables k1,k2,...kn from the MultiNest Hypercube to the physical space
        ! This is just a uniform prior, but has been sorted so that k1<k2<...<kn
        !
        if(Sorted%num_var /= 0 ) then
            temp_vec(1:Sorted%num_var) = Cube(Sorted%cube_pars)
            call  UniformSortPrior( temp_vec(1:Sorted%num_var), sorted_uniform_min,sorted_uniform_max)

            P(Sorted%cosmo_pars) = temp_vec(1:Sorted%num_var) 
        endif



    end subroutine Transform_Cube_to_Cosmo

    ! dumper, called after every updInt*10 iterations

    subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSlogZ, logZerr, context)
        !This routine is called after every updInt*10 iterations & at the end of the sampling allowing the posterior 
        !distribution & parameter constraints to be passed on to the user in the memory. The argument are as follows:
        !
        !The 2D arrays are Fortran arrays which are different to C/C++ arrays. In the example dumper routine provided with C & C++
        !eggbox examples, the Fortran arrays are copied on to C/C++ arrays.

        implicit none

        integer nSamples                ! number of samples in posterior array
        integer nlive                   ! number of live points
        integer nPar                    ! number of parameters saved (physical plus derived)
        double precision, pointer :: physLive(:,:) 
        !2D array containing the last set of live points (physical parameters
        !plus derived parameters) along with their loglikelihood values
        double precision, pointer :: posterior(:,:)
        !posterior distribution containing nSamples points. Each sample has nPar
        !parameters (physical + derived) along with the their loglike value &
        !posterior probability
        double precision, pointer :: paramConstr(:) ! array with mean, sigmas, maxlike & MAP parameters
        !     paramConstr(1, 1) to paramConstr(1, nPar)	     	 = mean values of the parameters
        !     paramConstr(1, nPar+1) to paramConstr(1, 2*nPar)   = standard deviation of the parameters
        !     paramConstr(1, nPar*2+1) to paramConstr(1, 3*nPar) = best-fit (maxlike) parameters
        !     paramConstr(1, nPar*4+1) to paramConstr(1, 4*nPar) = MAP (maximum-a-posteriori) parameters
        double precision maxLogLike         ! max loglikelihood value
        double precision logZ               ! log evidence
        double precision INSlogZ            ! log evidence as computed by Importance sampling 
        double precision logZerr            ! error on log evidence
        integer context                     ! not required by MultiNest, any additional information user wants to pass

        ! Note that at the moment this routine is not used in CosmoMC multinest

    end subroutine dumper

    !-----------------------------------------------------------------------

    subroutine nest_Sample

        implicit none

        integer context !total number of clusters found

        !Calling MultiNest
        call nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar, &
            nest_nCdims,maxClst,nest_updInt,nest_Ztol,nest_root,seed,nest_pWrap, &
            nest_fb,nest_resume,nest_outfile,initMPI,nest_logZero,nest_maxIter,getLogLikeNS,dumper,context)


    end subroutine nest_Sample

    !-------------------------------------------------------------------------




    !*************** Type(Prior) member functions ****************************

    !-------------------------------------------------------------------------
    subroutine initialize_prior(this,indx,list)
        ! This initialises all of the variables from an index defining the start
        ! of the region of the Cube where they will be stored and the list of
        ! indices to store in the cosmo_pars
        implicit none
        Class (Prior) :: this
        integer indx
        integer list(:)
        integer ii


        ! Assign the index
        this%cube_pos = indx

        ! Get the size of the array of parameter indices
        this%num_var = size(list)

        ! Allocate and define the cosmological parameters
        allocate(this%cosmo_pars(this%num_var))
        this%cosmo_pars= list

        ! Allocate and define the cube parameters
        allocate(this%cube_pars(this%num_var))
        do ii=1,this%num_var
            this%cube_pars(ii) = indx+ii
        enddo

    end subroutine initialize_prior

    subroutine add_point(this,val)
        ! This adds a point to cosmo_pars
        implicit none
        Class (Prior) :: this
        integer val
        integer temp_vec(max_num_params)
        integer temp_size

        if( allocated(this%cosmo_pars) ) then
            ! save current array
            temp_size=size(this%cosmo_pars)
            temp_vec(1:temp_size)=this%cosmo_pars

            ! resize cosmo_pars
            deallocate(this%cosmo_pars)
            allocate( this%cosmo_pars(temp_size+1) )

            ! restore saved values
            this%cosmo_pars(1:temp_size) = temp_vec(1:temp_size)

            ! add in the new value
            this%cosmo_pars(temp_size+1) = val

        else 
            ! If we need to allocate the array
            allocate(this%cosmo_pars(1))
            this%cosmo_pars(1) = val
        endif

        this%num_var = size(this%cosmo_pars)

    end subroutine add_point

    subroutine set_cube_pars(this,indx)
        ! This sets up the cube_pars vector
        implicit none
        Class (Prior) :: this
        integer indx
        integer ii

        ! set the index
        this%cube_pos = indx

        ! Erase the cube params if they're already allocated
        if( allocated(this%cube_pars) ) deallocate(this%cube_pars)

        ! Allocate and define the cube parameters
        ! Start with indx+1 and move to indx+num_var
        allocate(this%cube_pars(this%num_var))
        do ii=1,this%num_var
            this%cube_pars(ii) = indx+ii
        enddo

    end subroutine set_cube_pars




    subroutine setup_Gaussian_priors(prior_boundary,prior_vol)
        implicit none
        integer prior_boundary
        double precision prior_vol

        ! CosmoMC now has the ability to do gaussian priors
        ! The manner in which it does this is to add a compensating factor to
        ! the LogLikelihood.
        ! To be more in keeping with MultiNest, we should turn off this behavior
        ! by setting std=0, but record the values in a new variable

        allocate(nest_GaussPriors_std(size(BaseParams%GaussPriors%std)))
        nest_GaussPriors_std = BaseParams%GaussPriors%std

        allocate(nest_GaussPriors_mean(size(BaseParams%GaussPriors%mean)))
        nest_GaussPriors_mean = BaseParams%GaussPriors%mean

        BaseParams%GaussPriors%std = 0


    end subroutine setup_Gaussian_priors


    subroutine setup_MvG_priors(prior_boundary,prior_vol)
        implicit none
        integer prior_boundary
        double precision prior_vol

        integer temp_number                 ! number of parameters
        integer temp_params(max_num_params) ! vector of parameter indices
        double precision, dimension(:,:), allocatable :: pmat
        double precision det_covmat

        ! Get the indices of the named variables contained on the line
        ! 'MvG_parameters' and store them in temp_params
        ! also store the number of indices in temp_number (found by setting
        ! temp_number=-1)
        temp_number = -1
        call BaseParams%NameMapping%ReadIndices(trim(MvG_parameters), temp_params, temp_number)

        ! Set up the prior information using temp_params and temp_number
        ! found above
        call MvG%initialize_prior( prior_boundary, temp_params(1:temp_number) )
        ! Move the prior boundary on
        prior_boundary = prior_boundary + MvG%num_var

        ! Note, we will store the mean in the last column, and the
        ! covariance matrix in the first few
        allocate(MvG_prior_matrix(MvG%num_var+1,MvG%num_var))



        ! Read in the covariance matrix and its parameters
        allocate(pmat(max_num_params,max_num_params))
        call IO_ReadProposeMatrix(BaseParams%NameMapping,pmat,trim(file_MvG_cov))
        MvG_prior_matrix(1:MvG%num_var,1:MvG%num_var ) = pmat(MvG%cosmo_pars,MvG%cosmo_pars)

        ! Re-scale the covariance matrix by MvG_stretch
        MvG_prior_matrix(1:MvG%num_var,1:MvG%num_var)                     &
            =                                                               &
            MvG_stretch * MvG_prior_matrix(1:MvG%num_var,1:MvG%num_var)


        ! Read in the means
        call Read_Means( trim(file_MvG_mean), pmat(1,:) )
        MvG_prior_matrix(MvG%num_var+1,:) = pmat(1,MvG%cosmo_pars)


        ! Add this to the prior volume
        det_covmat = exp( MatrixSym_LogDet(MvG_prior_matrix(1:MvG%num_var,1:MvG%num_var)))
        prior_vol = prior_vol                      *   &
            pi**(MvG%num_var/2.d0)         /   &! volume of the...
            gamma(MvG%num_var/2.d0 + 1.d0) *   &! ...unit n-ball
            std_fac**(MvG%num_var)         *   &! expand volume to include 99% of probablity
            sqrt(det_covmat)


    end subroutine setup_MvG_priors

    subroutine setup_Ellipsoid_priors(prior_boundary,prior_vol)
        implicit none
        integer prior_boundary
        double precision prior_vol
        double precision det_covmat

        integer temp_number                 ! number of parameters
        integer temp_params(max_num_params) ! vector of parameter indices
        double precision, dimension(:,:), allocatable :: pmat

        ! Get the indices of the named variables contained on the line
        ! 'Ellipsoid_parameters' and store them in temp_params
        ! also store the number of indices in temp_number (found by setting
        ! temp_number=-1)
        temp_number = -1
        call BaseParams%NameMapping%ReadIndices(trim(Ellipsoid_parameters), temp_params, temp_number)

        ! Set up the prior information using temp_params and temp_number
        ! found above
        call Ellipsoid%initialize_prior( prior_boundary, temp_params(1:temp_number) )
        ! Move the prior boundary on
        prior_boundary = prior_boundary + Ellipsoid%num_var

        ! Note, we will store the mean in the last column, and the
        ! covariance matrix in the first few
        allocate(Ellipsoid_prior_matrix(Ellipsoid%num_var+1,Ellipsoid%num_var))



        ! Read in the covariance matrix and its parameters
        allocate(pmat(max_num_params,max_num_params))
        call IO_ReadProposeMatrix(BaseParams%NameMapping,pmat,trim(file_Ellipsoid_cov))
        Ellipsoid_prior_matrix(1:Ellipsoid%num_var,1:Ellipsoid%num_var ) = pmat(Ellipsoid%cosmo_pars,Ellipsoid%cosmo_pars)

        ! Re-scale the covariance matrix by Ellipsoid_stretch
        Ellipsoid_prior_matrix(1:Ellipsoid%num_var,1:Ellipsoid%num_var)     &
            =                                                               &
            Ellipsoid_stretch * Ellipsoid_prior_matrix(1:Ellipsoid%num_var,1:Ellipsoid%num_var)


        ! Read in the means
        call Read_Means( trim(file_Ellipsoid_mean), pmat(1,:) )
        Ellipsoid_prior_matrix(Ellipsoid%num_var+1,:) = pmat(1,Ellipsoid%cosmo_pars)


        ! Add this to the prior volume
        det_covmat = exp( MatrixSym_LogDet(Ellipsoid_prior_matrix(1:Ellipsoid%num_var,1:Ellipsoid%num_var)))
        prior_vol = prior_vol                    *   &
            pi**(Ellipsoid%num_var/2.d0)         /   &! volume of the...
            gamma(Ellipsoid%num_var/2.d0 + 1.d0) *   &! ...unit n-ball
            sqrt(det_covmat)


    end subroutine setup_Ellipsoid_priors

    subroutine setup_Sorted_priors(prior_boundary,prior_vol)
        implicit none
        integer prior_boundary
        double precision prior_vol

        integer temp_number                 ! number of parameters
        integer temp_params(max_num_params) ! vector of parameter indices



        ! Get the indices of the named variables contained on the line
        ! 'sorted_uniform_parameters' and store them in temp_params
        ! also store the number of indices in temp_number (found by setting
        ! temp_number=-1)
        temp_number = -1
        call BaseParams%NameMapping%ReadIndices(trim(sorted_uniform_parameters), temp_params, temp_number)

        ! Set up the prior information using temp_params and temp_number
        ! found above
        call Sorted%initialize_prior( prior_boundary, temp_params(1:temp_number) )
        ! Move the prior boundary on
        prior_boundary = prior_boundary + Sorted%num_var

        ! The prior volume increase is just by a n dimensional triangle
        prior_vol = prior_vol * (sorted_uniform_max - sorted_uniform_min)**Sorted%num_var/ gamma(1.d0 + Sorted%num_var)

    end subroutine setup_Sorted_priors

    subroutine setup_Separable_priors(prior_boundary,prior_vol)
        implicit none
        integer prior_boundary
        double precision prior_vol
        integer i
        logical i_assigned ! switch to tell if we've assigned it a type

        do i=1,num_params_used
            i_assigned=.false.

            if ( MvG%num_var/=0) then
                    ! If it's a Multivariate gaussian prior, then skip it
                if ( any( params_used(i) == MvG%cosmo_pars ) ) i_assigned=.true.
            endif

            if ( Ellipsoid%num_var/=0 ) then
                    ! If it's a Ellipsoid uniform prior, then skip it
                if (any( params_used(i) == Ellipsoid%cosmo_pars ) )  i_assigned=.true.
            endif

            if ( Sorted%num_var/=0 ) then
                ! If it's a Sorted uniform prior, then skip it
                if (any( params_used(i) == Sorted%cosmo_pars ) ) i_assigned=.true.
            endif

            if (0d0 /= nest_GaussPriors_std(params_used(i)) ) then 
                ! If its a standard cosmomc Gaussian prior paramter
                ! Add this integer to the Gaussian cosmo_pars array
                call Gaussian%add_point( params_used(i) )

                ! Add this to the prior volume
                prior_vol = prior_vol * std_fac *  2*nest_GaussPriors_std(params_used(i))
                i_assigned=.true.
            endif

            if (i_assigned==.false.) then
                ! Otherwise, assume uniform
                ! Add this integer to the Uniform cosmo_pars array
                call Uniform%add_point( params_used(i) )

                ! Add this to the prior volume
                prior_vol = prior_vol * &
                    (  BaseParams%Pmax(params_used(i)) - BaseParams%Pmin(params_used(i))  )
            endif

        enddo

        ! Set up the cube_pars for use later
        call Gaussian%set_cube_pars(prior_boundary)
        prior_boundary = prior_boundary + Gaussian%num_var
        call Uniform%set_cube_pars(prior_boundary)
        prior_boundary = prior_boundary + Uniform%num_var


    end subroutine setup_Separable_priors


    subroutine Read_Means( file_mean, out_vec )
        implicit none

        character(LEN=*), intent(in) :: file_mean
        double precision, intent(out) :: out_vec(:)
        character(LEN=ParamNames_maxlen) name
        double precision mean
        integer ix, status
        Type(TTextFile) :: F

        call F%Open(file_mean)
        do
            read(F%unit, *, iostat=status) name, mean
            ix = BaseParams%NameMapping%index(name)
            if (ix/=-1) out_vec(ix) = mean
            if (status/=0) exit
        end do
        call F%Close()

    end subroutine Read_Means

    subroutine Trim_Volume(ndims,nPar,prior_vol)
        use RandUtils
        implicit none
        integer, intent(in) :: ndims,nPar
        double precision :: prior_vol
        ! The fact that multinest doesn't use points outside physical ranges
        ! means that the prior volume needs to be trimmed down by some factor.
        !
        ! This function trims the volume by calculating that factor.
        ! It draws random variables inside the prior range, and records which
        ! ones are 'in' the physical range.
        !
        ! We then estimate the factor and its error using a bayesian estimator
        ! based on a logarithmic prior on the factor:
        !
        ! Let f be the factor of the volume lying within the physically valid
        ! prior range. We draw some samples, recording num_in inside the region
        ! and num_out outside the region. P(num_in,num_out|f) is therefore
        ! binomial. With a logarithmic prior, the posterior P(f|num_in,num_out) 
        ! is a distribution with mean and variance:
        !
        ! <f>    = num_in / ( 1 + num_in + num_out )
        ! var(f) = num_in * ( 1 + num_out) / 
        !          ( (1+num_in+num_out)^2(2+num_in+num_out) )
        !
        ! We draw samples until sqrt(var(f))/<f> < 0.01, i.e. the fractional
        ! error in our estimate of f is less than 1%.
        
        integer, parameter ::  reduce_cycle = 1000

        integer j
        integer num_in_p, num_out_p
        integer num_in, num_out

        double precision Cube(nPar), P(nPar)
        logical still_calculating
#ifdef MPI
        integer ierror
#endif

        num_out_p=0
        num_in_p=0
        still_calculating=.true.


        do while (still_calculating)

            ! Generate a point in the cube
            do j=1,ndims
                Cube(j)=ranmar()
            enddo

            ! Transform to the cosmological parameter space
            call Transform_Cube_to_Cosmo(Cube,P,nPar)

            ! Record whether this outside or inside the physical parameter space
            if (any( P(params_used(:num_params_used)) < BaseParams%Pmin(params_used(:num_params_used))) .or. any( P(params_used(:num_params_used)) > BaseParams%Pmax(params_used(:num_params_used))) ) then
                num_out_p=num_out_p+1
            else
                num_in_p=num_in_p+1
            endif

            ! We reduce when each core has done 1000 iterations 
            if ( mod(num_in_p + num_out_p,reduce_cycle) == 0 ) then
                ! Find the total number in and out, by summing across nodes and broadcast it to all nodes
                call MPI_Allreduce( num_in_p , num_in , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror )
                call MPI_Allreduce( num_out_p, num_out, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror )

                ! Test to see whether we're within the error margin
                if ( dble(1+num_out) < prior_vol_frac_err*prior_vol_frac_err *(dble(num_in)*dble(2+num_out+num_in)) ) still_calculating = .false.
                
            endif

        enddo

        ! trim the volume down by the fraction that are inside the physical
        ! parameter space
        prior_vol = dble(num_in)/dble(1+num_in+num_out) * prior_vol

        if ( IsMainMPI() ) &
            write(*,'("Prior Volume needed to be trimmed by a factor of :",1e10.3e2,"(1+/-",1f5.2,")")') &
            dble(num_in)/dble(1+num_in+num_out),prior_vol_frac_err


    end subroutine Trim_Volume














end module nestwrap
