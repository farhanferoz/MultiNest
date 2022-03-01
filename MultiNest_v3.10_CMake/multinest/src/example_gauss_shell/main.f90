program main

	use params
      	use nestwrapper
      
	implicit none
      
      
	integer i
      
	!setting priors
	spriorran(1:sdim,1)=-6.
      	spriorran(1:sdim,2)=6.
      
      	!no parameters to wrap around
      	nest_pWrap=0
      
      	!setting the centers for dimensions other than the first one
      	do i=1,sModes
      		sc(i,2:sdim)=0.
      	end do

      	call nest_Sample
      	stop
end
