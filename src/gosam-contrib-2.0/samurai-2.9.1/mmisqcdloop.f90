module mmisqcdloop
       use precision, only: ki, ki_ql
       use constants
       use options
       use mfunctions
       use notfirst
       use mmishighrank
       implicit none
    
       private
    
       complex(ki), dimension(:,:), allocatable, public :: s_mat
    
    	interface
    		function qlI4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
    			use precision, only: ki_ql
    			implicit none
    			real(ki_ql), intent(in) :: p1,p2,p3,p4,s12,s23
    			real(ki_ql), intent(in) :: m1,m2,m3,m4,mu2
    			integer, intent(in) :: ep
    			complex(ki_ql) :: qlI4
    		end function qlI4
    	end interface
    	interface
    		function qlI3(p1,p2,p3,m1,m2,m3,mu2,ep)
    			use precision, only: ki_ql
    			implicit none
    			real(ki_ql), intent(in) :: p1,p2,p3
    			real(ki_ql), intent(in) :: m1,m2,m3,mu2
    			integer, intent(in) :: ep
    			complex(ki_ql) :: qlI3
    		end function qlI3
    	end interface
    	interface
    		function qlI2(p1,m1,m2,mu2,ep)
    			use precision, only: ki_ql
    			implicit none
    			real(ki_ql), intent(in) :: p1
    			real(ki_ql), intent(in) :: m1,m2,mu2
    			integer, intent(in) :: ep
    			complex(ki_ql) :: qlI2
    		end function qlI2
    	end interface
    	interface
    		function qlI1(m1,mu2,ep)
    			use precision, only: ki_ql
    			implicit none
    			real(ki_ql), intent(in) :: m1,mu2
    			integer, intent(in) :: ep
    			complex(ki_ql) :: qlI1
    		end function qlI1
    	end interface
    
    	public :: qcdloop4, qcdloop3, qcdloop2,qcdloop2hr, qcdloop1
    
    contains
    
    subroutine qcdloop4(V,m,scale2,MI4,cache_flag, cache_offset, scalar_cache)
    	implicit none
    	real(ki),    dimension(1:6),  intent(in ) :: V
    	complex(ki), dimension(0:3),  intent(in ) :: m
    	real(ki), 		      intent(in ) :: scale2
    	complex(ki), dimension(-2:0), intent(out) :: MI4
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1, m2, m3
    	real(ki)    :: V1, V2, V3, V21, V31, V32
    	integer     :: ep, cache_index
    	complex(ki) :: MI4tmp
    	
    	notfirsti=.true.
    	
    	m0  = m(0)
    	m1  = m(1)
    	m2  = m(2)
    	m3  = m(3)
    	V1  = V(1)
    	V2  = V(2)
    	V3  = V(3)
    	V21 = V(4)
    	V31 = V(5)
    	V32 = V(6)
    	
      1   Format(A3,I4,A1,I2,A5,D24.15,A1,D24.15,A3)
    
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	if ( (aimag(m0).ne.0) .or. (aimag(m1).ne.0) .or. &
    	&    (aimag(m2).ne.0) .or. (aimag(m3).ne.0) ) then
    		print *,'isca=1: QCDLoop does not support complex masses'
    		stop
    	endif
    	do ep=-2,0
    		if (present(cache_flag)) then
    			if (cache_flag) then
    				MI4tmp = scalar_cache(ep,cache_index)
    			else
    				MI4tmp=qlI4(&
    				& real(V1,ki_ql),real(V21,ki_ql),real(V32,ki_ql),&
    				& real(V3,ki_ql),real(V2,ki_ql),real(V31,ki_ql),&
    				& real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
    				& real(m3,ki_ql),real(scale2,ki_ql),ep)
    				scalar_cache(ep,cache_index) = MI4tmp
    			end if
    		else
    			MI4tmp=qlI4(&
    			& real(V1,ki_ql),real(V21,ki_ql),real(V32,ki_ql),&
    			& real(V3,ki_ql),real(V2,ki_ql),real(V31,ki_ql),&
    			& real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
    			& real(m3,ki_ql),real(scale2,ki_ql),ep)
    		end if
    		MI4(ep)=MI4tmp
    	enddo
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    end subroutine qcdloop4
    
    subroutine qcdloop3(V,m,scale2,MI3,cache_flag, cache_offset, scalar_cache)
    	implicit none
    	
    	real(ki),    dimension(1:3),  intent(in ) :: V
    	complex(ki), dimension(0:2),  intent(in ) :: m
    	real(ki), 		      intent(in ) :: scale2
    	complex(ki), dimension(-2:0), intent(out) :: MI3
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1, m2
    	real(ki)    :: V1, V2, V3
    	complex(ki) :: MI3tmp
    	integer     :: ep, cache_index
    	
    	notfirsti=.true.
    
    	m0 = m(0)
    	m1 = m(1)
    	m2 = m(2)
    	V1 = V(1)
    	V2 = V(2)
    	V3 = V(3)
    	
      1  Format(A3,I3,A1,I2,A5,D24.15,A1,D24.15,A3)
    
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	do ep=-2,0
    		if (present(cache_flag)) then
    			if (cache_flag) then
    				MI3tmp = scalar_cache(ep,cache_index)
    			else
    				MI3tmp=qlI3(&
    				& real(V1,ki_ql),real(V2,ki_ql),real(V3,ki_ql),&
    				& real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
    				& real(scale2,ki_ql),ep)
    				scalar_cache(ep,cache_index) = MI3tmp 
    			end if
    		else
    			MI3tmp=qlI3(&
    			& real(V1,ki_ql),real(V2,ki_ql),real(V3,ki_ql),&
    			& real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
    			& real(scale2,ki_ql),ep)
    		end if
    		MI3(ep)=MI3tmp
    	enddo
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    
    end subroutine qcdloop3
    
    subroutine qcdloop2(K11,m,scale2,MI2J0,MI2J1,MI2J00,MI2J01,MI2J11,&
    			& cache_flag, cache_offset, scalar_cache)
    	implicit none
    	
    	real(ki),    		     intent(in ) :: K11 
    	complex(ki), dimension(0:1), intent(in ) :: m
    	real(ki), 		     intent(in ) :: scale2
    	complex(ki), dimension(-2:0),intent(out) :: MI2J0, MI2J1, MI2J00, MI2J01, MI2J11
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1
    	integer     :: ep, cache_index
    
    	notfirsti=.true.
    	
    	m0=m(0)
    	m1=m(1)
    	
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	if (present(cache_flag)) then
    		if (cache_flag) then
    			MI2J0(:)  = scalar_cache(:,cache_index+0)
    			MI2J1(:)  = scalar_cache(:,cache_index+1)
    			MI2J01(:) = scalar_cache(:,cache_index+2)
    			MI2J11(:) = scalar_cache(:,cache_index+3)
    			MI2J00(:) = scalar_cache(:,cache_index+4)
    		else
    			 do ep=-2,0
    				MI2J00(ep) = qlI2(&
    				& real(K11,ki_ql),real(m0,ki_ql),real(m0,ki_ql),&
    				& real(scale2,ki_ql),ep)
    				MI2J11(ep) = qlI2(&
    				& real(K11,ki_ql),real(m1,ki_ql),real(m1,ki_ql),&
    				& real(scale2,ki_ql),ep)
    				MI2J01(ep)=qlI2(&
    				& real(K11,ki_ql),real(m0,ki_ql),real(m1,ki_ql),&
    				& real(scale2,ki_ql),ep)
    				MI2J0(ep) = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
    				MI2J1(ep) = qlI1(real(m1,ki_ql),real(scale2,ki_ql),ep)
    			 end do
    			 scalar_cache(:,cache_index+0) = MI2J0(:)
    			 scalar_cache(:,cache_index+1) = MI2J1(:)
    			 scalar_cache(:,cache_index+2) = MI2J01(:)
    			 scalar_cache(:,cache_index+3) = MI2J11(:)
    			 scalar_cache(:,cache_index+4) = MI2J00(:)
    		end if
    	else
    		do ep=-2,0
    			MI2J00(ep) = qlI2(&
    			& real(K11,ki_ql),real(m0,ki_ql),real(m0,ki_ql),&
    			& real(scale2,ki_ql),ep)
    			MI2J11(ep) = qlI2(&
    			& real(K11,ki_ql),real(m1,ki_ql),real(m1,ki_ql),&
    			& real(scale2,ki_ql),ep)
    			MI2J01(ep)=qlI2(&
    			& real(K11,ki_ql),real(m0,ki_ql),real(m1,ki_ql),&
    			& real(scale2,ki_ql),ep)
    			MI2J0(ep) = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
    			MI2J1(ep) = qlI1(real(m1,ki_ql),real(scale2,ki_ql),ep)
    		end do
    	end if
    	if (present(cache_flag)) cache_offset = cache_offset + 5
    
    end subroutine qcdloop2
    
    
    subroutine qcdloop2hr(K11,m,scale2,J111)
    	implicit none
    	real(ki),    intent(in) :: K11
    	complex(ki), dimension(0:1), intent(in) :: m
    	real(ki),    intent(in) :: scale2
    	complex(ki), dimension(-2:0), intent(out) :: J111
    	integer     :: ep
    	complex(ki) :: m0, m1
     	complex(ki), dimension(-2:0) :: B0p12, B0z11, B0z22
    
    	m0 = m(0)
    	m1 = m(1)
    
    	do ep=-2,0
    		B0p12(ep) = qlI2(&
    		& real(K11,ki_ql),real(m0,ki_ql),real(m1,ki_ql),&
    		& real(scale2,ki_ql),ep)
    		B0z11(ep) = qlI2(&
    		& real(zip,ki_ql),real(m0,ki_ql),real(m0,ki_ql),&
    		& real(scale2,ki_ql),ep)
    		B0z22(ep) = qlI2(&
    		& real(zip,ki_ql),real(m1,ki_ql),real(m1,ki_ql),&
    		& real(scale2,ki_ql),ep)
    	enddo
    
    	call HJ111(J111, K11,m0,m1,B0p12,B0z11,B0z22)
    
    end subroutine qcdloop2hr
    
    subroutine qcdloop1(m,scale2,MI1, cache_flag, cache_offset, scalar_cache)
    	implicit none 
    	complex(ki), 			  intent(in ) :: m
    	real(ki), 			  intent(in ) :: scale2
    	complex(ki), dimension(-2:0), 	  intent(out) :: MI1
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	integer     :: j1
    	complex(ki) :: m0, MI1tmp
    	integer     :: ep, cache_index
    	
    	notfirsti=.true.
    	
    	m0=m
    	
      1   Format(A3,I2,A1,I2,A5,D24.15,A1,D24.15,A3)
    
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	do ep=-2,0
    		if (present(cache_flag)) then
    			if (cache_flag) then
    				MI1tmp = scalar_cache(ep,cache_index)
    			else
    				MI1tmp = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
    				scalar_cache(ep,cache_index) = MI1tmp
    			end if
    		else
    			MI1tmp = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
    		end if
    		MI1(ep)=MI1tmp
    	enddo
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    
    end subroutine qcdloop1
    
end module mmisqcdloop

