module mmisgolem 
    	use precision, only: ki
    	use precision_golem, only: ki_gol => ki
    	use constants
    	use options
    	use notfirst
    	use mmishighrank
    	implicit none
    	
    	private
    	
    	interface
    		function gD0(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			real(ki), intent(in) :: p1,p2,p3,p4,s12,s23
    			real(ki), intent(in) :: m1,m2,m3,m4
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gD0
    		end function gD0
    	end interface
    	interface
    		function gC0(p1,p2,p3,m1,m2,m3,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			real(ki), intent(in) :: p1,p2,p3
    			real(ki), intent(in) :: m1,m2,m3
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gC0
    		end function gC0
    	end interface
    	interface
    		function gB0(p1,m1,m2,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			real(ki), intent(in) :: p1
    			real(ki), intent(in) :: m1,m2
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gB0
    		end function gB0
    	end interface
    	interface
    		function gD0C(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			complex(ki), intent(in) :: p1,p2,p3,p4,s12,s23
    			complex(ki), intent(in) :: m1,m2,m3,m4
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gD0C
    		end function gD0C
    	end interface
    	interface
    		function gC0C(p1,p2,p3,m1,m2,m3,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			complex(ki), intent(in) :: p1,p2,p3
    			complex(ki), intent(in) :: m1,m2,m3
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gC0C
    		end function gC0C
    	end interface
    	interface
    		function gC0i(idt,p1,p2,p3,m1,m2,m3,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			real(ki), intent(in) :: p1,p2,p3
    			real(ki), intent(in) :: m1,m2,m3
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gC0i
    		end function gC0i
    	end interface
    	interface
    		function gC0iC(idt,p1,p2,p3,m1,m2,m3,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			complex(ki), intent(in) :: p1,p2,p3
    			complex(ki), intent(in) :: m1,m2,m3
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gC0iC
    		end function gC0iC
    	end interface
    	interface
    		function gB0C(p1,m1,m2,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			complex(ki), intent(in) :: p1
    			complex(ki), intent(in) :: m1,m2
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gB0C
    		end function gB0C
    	end interface
    	interface
    		function gB0i(idt,p1,m1,m2,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			real(ki), intent(in) :: p1
    			real(ki), intent(in) :: m1,m2
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gB0i
    		end function gB0i
    	end interface
    	interface
    		function gB0iC(idt,p1,m1,m2,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			complex(ki), intent(in) :: p1
    			complex(ki), intent(in) :: m1,m2
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gB0iC
    		end function gB0iC
    	end interface
    	!interface
    	!	function gA0(m1,mu2,ep)
    	!		use precision_golem, only: ki
    	!		implicit none
    	!		real(ki), intent(in) :: m1
    	!		real(ki), intent(in) :: mu2
    	!		integer, intent(in) :: ep
    	!		complex(ki) :: gA0
    	!	end function gA0
    	!end interface
    	interface
    		function gA0C(m1,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			complex(ki), intent(in) :: m1
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gA0C
    		end function gA0C
    	end interface
    	interface
    		function gA0i(idt,m1,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			real(ki), intent(in) :: m1
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gA0i
    		end function gA0i
    	end interface
    	interface
    		function gA0iC(idt,m1,mu2,ep)
    			use precision_golem, only: ki
    			implicit none
    			character (len=*), intent (in) :: idt
    			complex(ki), intent(in) :: m1
    			real(ki), intent(in) :: mu2
    			integer, intent(in) :: ep
    			complex(ki) :: gA0iC
    		end function gA0iC
    	end interface
    
    
    	public :: golemMI4, golemMI3, golemMI2, golemMI2hr1,golemMI2hr2, golemMI1
    
    
    contains
    
    subroutine golemMI4(V,m,scale2,MI4,cache_flag, cache_offset, scalar_cache)
    	implicit none
    	real(ki),    dimension(1:6), intent(in ) :: V
    	complex(ki), dimension(0:3), intent(in ) :: m
    	real(ki),		     intent(in ) :: scale2
    	complex(ki), dimension(-2:0),       intent(out) :: MI4
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1, m2, m3
    	complex(ki) :: V1, V2, V3, V21, V31, V32
    	integer     :: ep, cache_index
    
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
    
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    
    	call gtrunc(abs(V32)+abs(V31),V1,V2,V3,V21,V32,V31,m0,m1,m2,m3)
    	do ep=-2,0
    		if (present(cache_flag)) then
    			if (cache_flag) then
    				MI4(ep) = scalar_cache(ep,cache_index)
    			else
    				MI4(ep)=gD0(real(V1,ki_gol),real(V21,ki_gol),&
    				& real(V32,ki_gol),real(V3,ki_gol),&
    				& real(V2,ki_gol),real(V31,ki_gol),&
    				& real(m0,ki_gol),real(m1,ki_gol),&
    				& real(m2,ki_gol),real(m3,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				scalar_cache(ep,cache_index) = MI4(ep)
    			end if
    		else
    			MI4(ep)=gD0(real(V1,ki_gol),real(V21,ki_gol),&
    			& real(V32,ki_gol),real(V3,ki_gol),&
    			& real(V2,ki_gol),real(V31,ki_gol),&
    			& real(m0,ki_gol),real(m1,ki_gol),&
    			& real(m2,ki_gol),real(m3,ki_gol),&
    			& real(scale2,ki_gol),ep)
    		end if
    	end do
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    end subroutine golemMI4
    
    subroutine golemMI3(V,m,scale2,MI3,cache_flag, cache_offset, scalar_cache)
    	implicit none
    	real(ki),    dimension(1:3),  intent(in ) :: V
    	complex(ki), dimension(0:2),  intent(in ) :: m
    	real(ki), 		      intent(in ) :: scale2
    	complex(ki), dimension(-2:0), intent(out) :: MI3
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1, m2
    	complex(ki) :: V1, V2, V3
    	integer :: ep, cache_index
    	complex(ki), dimension(-2:0) :: c0t
    	
    	notfirsti=.true.
    
    	m0 = m(0)
    	m1 = m(1)
    	m2 = m(2)
    	V1 = V(1)
    	V2 = V(2)
    	V3 = V(3)
    	
    
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	call gtrunc(abs(V1)+abs(V2)+abs(V3),V1,V2,V3,m0,m1,m2)
    	do ep=-2,0
    		if (present(cache_flag)) then
    			if (cache_flag) then
    				c0t(ep) = scalar_cache(ep,cache_index)
    			else
    				c0t(ep)=gC0(real(V1,ki_gol),real(V2,ki_gol),&
    				& real(V3,ki_gol),real(m0,ki_gol),&
    				& real(m1,ki_gol),real(m2,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				scalar_cache(ep,cache_index) = c0t(ep)
    			end if
    		else
    			c0t(ep)=gC0(real(V1,ki_gol),real(V2,ki_gol),&
    			& real(V3,ki_gol),real(m0,ki_gol),&
    			& real(m1,ki_gol),real(m2,ki_gol),&
    			& real(scale2,ki_gol),ep)
    		end if
    	end do
    	MI3(:) = c0t(:) 
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    end subroutine golemMI3
    
    subroutine golemMI2(K11in,m,scale2,J0,J1,J00,J01,J11, cache_flag, cache_offset, scalar_cache)
    	implicit none
    	real(ki),     		      intent(in ) :: K11in
    	complex(ki), dimension(0:1),  intent(in ) :: m
    	real(ki), 		      intent(in ) :: scale2
            complex(ki), dimension(-2:0),intent(out) :: J0, J1, J00, J01, J11
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	complex(ki) :: m0, m1
    	complex(ki) :: K11
    	integer     :: ep, cache_index
    
    	notfirsti=.true.
    	
    	m0=m(0)
    	m1=m(1)
    	K11=K11in
    	
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	if (present(cache_flag)) then
    		if (cache_flag) then
    			J0(:)  = scalar_cache(:,cache_index+0)
    			J1(:)  = scalar_cache(:,cache_index+1)
    			J01(:) = scalar_cache(:,cache_index+2)
    			J11(:) = scalar_cache(:,cache_index+3)
    			J00(:) = scalar_cache(:,cache_index+4)
    		else
    			call gtrunc(abs(K11)+1.0_ki,K11,m0,m1)
    			do ep=-2,0
    				J00(ep)= gB0(real(K11,ki_gol),&
    				& real(m0,ki_gol),real(m0,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				J11(ep)= gB0(real(K11,ki_gol),&
    				& real(m1,ki_gol),real(m1,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				J01(ep)= gB0(real(K11,ki_gol),&
    				& real(m0,ki_gol),real(m1,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				J0(ep) = gA0(real(m0,ki_gol),&
    				& real(scale2,ki_gol),ep)
    				J1(ep) = gA0(real(m1,ki_gol),&
    				& real(scale2,ki_gol),ep)
    			end do
    			scalar_cache(:,cache_index+0) = J0(:)
    			scalar_cache(:,cache_index+1) = J1(:)
    			scalar_cache(:,cache_index+2) = J01(:)
    			scalar_cache(:,cache_index+3) = J11(:)
    			scalar_cache(:,cache_index+4) = J00(:)
    		end if
    	else
    		call gtrunc(abs(K11)+1.0_ki, K11,m0,m1)
    		do ep=-2,0
    			J00(ep)= gB0(real(K11,ki_gol),&
    			& real(m0,ki_gol),real(m0,ki_gol),&
    			& real(scale2,ki_gol),ep)
    			J11(ep)= gB0(real(K11,ki_gol),&
    			& real(m1,ki_gol),real(m1,ki_gol),&
    			& real(scale2,ki_gol),ep)
    			J01(ep)= gB0(real(K11,ki_gol),&
    			& real(m0,ki_gol),real(m1,ki_gol),&
    			& real(scale2,ki_gol),ep)
    			J0(ep) = gA0(real(m0,ki_gol),&
    			& real(scale2,ki_gol),ep)
    			J1(ep) = gA0(real(m1,ki_gol),&
    			& real(scale2,ki_gol),ep)
    		end do
    	end if
    	if (present(cache_flag)) cache_offset = cache_offset + 5
    end subroutine golemMI2
    
    subroutine golemMI2hr1(K11,m,scale2,J111,J001)
    	implicit none
    	real(ki),     		      intent(in )           :: K11
    	complex(ki), dimension(0:1),  intent(in )           :: m
    	real(ki), 		      intent(in )           :: scale2
            complex(ki), dimension(-2:0),intent(out)           :: J111
            complex(ki), dimension(-2:0),intent(out), optional :: J001
    	
    	complex(ki) :: m0, m1
    	integer     :: ep
    
    	m0=m(0)
    	m1=m(1)
    
    		do ep=-2,0
    			if(present(J001)) then
    				J001(ep) = gB0i("001",real(K11,ki_gol),&
    				& real(m0,ki_gol),real(m1,ki_gol),&
    				& real(scale2,ki_gol),ep)
    			end if
    			J111(ep) = gB0i("111",real(K11,ki_gol),&
    			& real(m0,ki_gol),real(m1,ki_gol),&
    			& real(scale2,ki_gol),ep)
    		end do
    end subroutine golemMI2hr1
    
    subroutine golemMI2hr2(K11,m,scale2,J111)
    	implicit none
    	real(ki),   		      intent(in ) :: K11
    	complex(ki), dimension(0:1),  intent(in ) :: m
    	real(ki),   		      intent(in ) :: scale2
    	complex(ki), dimension(-2:0), intent(out) :: J111
    
    	integer     :: ep
    	complex(ki) :: m0,m1
     	complex(ki), dimension(-2:0) :: B0p12, B0z11, B0z22
    
    	m0=m(0)
    	m1=m(1)
    
    	do ep=-2,0
    		B0p12(ep)= gB0(real(K11,ki_gol),&
    		& real(m0,ki_gol),real(m1,ki_gol),&
    		& real(scale2,ki_gol),ep)
    		B0z11(ep)= gB0(real(zip,ki_gol),&
    		& real(m0,ki_gol),real(m0,ki_gol),&
    		& real(scale2,ki_gol),ep)
    		B0z22(ep)= gB0(real(zip,ki_gol),&
    		& real(m1,ki_gol),real(m1,ki_gol),&
    		& real(scale2,ki_gol),ep)
    	enddo
    
    	call HJ111(J111, K11,m0,m1,B0p12,B0z11,B0z22)
    
    end subroutine golemMI2hr2
    
    
    
    subroutine golemMI1(m,scale2,MI1,cache_flag, cache_offset, scalar_cache)
    
    	implicit none 
    	complex(ki), 			  intent(in ) :: m
    	real(ki), 			  intent(in ) :: scale2
    	complex(ki), dimension(-2:0),     intent(out) :: MI1
    	
    	logical,     intent(in   ), optional 			:: cache_flag
    	integer,     intent(inout), optional 			:: cache_offset
    	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
    	
    	integer :: ep, cache_index
    	complex(ki) :: m0
    
    	m0=m
    	
    	notfirsti=.true.
    	
    	MI1(-2)=czip
    	if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
    	if (present(cache_flag)) then
    		if (cache_flag) then
    			do ep=-1,0
    				MI1(ep) = scalar_cache(ep,cache_index)
    			enddo
    		else
    			scalar_cache(-2,cache_index) = czip
    			call gtrunc(1.0_ki, m0)
    			do ep=-1,0
    				MI1(ep) = gA0(real(m0,ki_gol),real(scale2,ki_gol),ep)
    				scalar_cache(ep,cache_index) = MI1(ep)
    			enddo
    		end if
    	else
    		call gtrunc(1.0_ki, m0)
    		do ep=-1,0
    			MI1(ep) = gA0(real(m0,ki_gol),real(scale2,ki_gol),ep)
    		enddo
    	end if
    	if (present(cache_flag)) cache_offset = cache_offset + 1
    end subroutine golemMI1
    
    
    pure subroutine gtrunc(ref,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
    	implicit none
    	real(ki), intent(in) :: ref
    	complex(ki), intent(inout), optional :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
    	real(ki), parameter :: small = 1.0E-08_ki
    	if(present(s1)) then
    		if(abs(s1/ref) .lt. small) s1 = 0.0_ki
    	end if
    	if(present(s2)) then
    		if(abs(s2/ref) .lt. small) s2 = 0.0_ki
    	end if
    	if(present(s3)) then
    		if(abs(s3/ref) .lt. small) s3 = 0.0_ki
    	end if
    	if(present(s4)) then
    		if(abs(s4/ref) .lt. small) s4 = 0.0_ki
    	end if
    	if(present(s5)) then
    		if(abs(s5/ref) .lt. small) s5 = 0.0_ki
    	end if
    	if(present(s6)) then
    		if(abs(s6/ref) .lt. small) s6 = 0.0_ki
    	end if
    	if(present(s7)) then
    		if(abs(s7/ref) .lt. small) s7 = 0.0_ki
    	end if
    	if(present(s8)) then
    		if(abs(s8/ref) .lt. small) s8 = 0.0_ki
    	end if
    	if(present(s9)) then
    		if(abs(s9/ref) .lt. small) s9 = 0.0_ki
    	end if
    	if(present(s10)) then
    		if(abs(s10/ref) .lt. small) s10 = 0.0_ki
    	end if
    end  subroutine gtrunc
    
    function gA0(m0,mu2,ep)
       implicit none
       real(ki_gol), intent(in) :: m0, mu2
       integer, intent(in) :: ep
       complex(ki_gol) :: gA0
       if(ep.eq.(-2) .or. m0.eq.0.0_ki_gol) then
          gA0 = (0.0_ki_gol, 0.0_ki_gol)
       elseif(ep.eq.(-1)) then
          gA0 = m0 * gB0(0.0_ki_gol,m0,m0,mu2,-1)
       else
          gA0 = m0 * (gB0(0.0_ki_gol,m0,m0,mu2,0) &
              &    +  gB0(0.0_ki_gol,m0,m0,mu2,-1))
       end if
    end function gA0


end module mmisgolem

