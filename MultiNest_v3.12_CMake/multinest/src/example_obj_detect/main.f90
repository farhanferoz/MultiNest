program main

	use params
	use nestwrapper
      
      	implicit none
	integer i
      	double precision temp(4),s_nullev
      
      	sdim=0
	do i=1,4
      		if(spriorran(i,1)/=spriorran(i,2)) sdim=sdim+1
      	end do
      	nest_nPar=sdim
	nest_pWrap=0
      
  	call sCreateData(temp,sodata,.false.) 
      
 	call sInitialize(s_nullev)
      	call nest_Sample 
      
      	write(*,*)"Null Evidence=",s_nullev   
      
      	stop
end
