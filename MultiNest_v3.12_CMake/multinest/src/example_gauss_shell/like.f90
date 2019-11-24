module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	double precision temp(sdim),dist,loclik
	integer i,j
	double precision TwoPi
         
	TwoPi=6.2831853

	slhood=-huge(1.d0)*epsilon(1.d0)
	
	!rescaling the parameters in unit hypercube according to the prior
	do i=1,sdim
		temp(i)=(spriorran(i,2)-spriorran(i,1))*Cube(i)+spriorran(i,1)
	end do
	Cube(1:sdim)=temp(1:sdim)

	do i=1,sModes
		dist=(sqrt(sum((temp(1:sdim)-sc(i,1:sdim))**2.))-sr(i))**2
		loclik=-dist/(2.*(sw(i)**2.))-log(TwoPi*sw(i)**2)/2.
		slhood=logSumExp(slhood,loclik)
	end do

end subroutine slikelihood
      
!=======================================================================

end module like
