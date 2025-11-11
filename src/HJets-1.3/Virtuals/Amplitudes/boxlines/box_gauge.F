c###################################################################
      subroutine box_gaugeOneLoop(p1,p2,p3,p4,musq,result,resultb,ward)
C Compute the two gauge tests for the boxes and result. .True. if the max of it
c is smaller than 1d-2. 

      IMPLICIT NONE
      DOUBLE PRECISION p1(0:3), p2(0:3), p3(0:3), p4(0:3) 
      DOUBLE COMPLEX result(2), resultb(2)
      DOUBLE PRECISION musq, testb(2)
      DOUBLE COMPLEX cte
      DOUBLE PRECISION dotrr,pi,pi2
      EXTERNAL dotrr
      PARAMETER (pi =3.14159265358979323846264338327950288d0)
      PARAMETER (pi2 =9.86960440108935861883449099987615114d0)
      DOUBLE PRECISION s14,lm,ls,ths,theta
      LOGICAL ward
      double complex  lns
      external theta
      double precision Gaugethres,gaugeAccuracy(2)
      common/Gauge/gaugethres,gaugeAccuracy

      ward=.False.

      s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
      lns=log(Abs(s14))-(0,1d0)*pi*theta(s14)-log(musq)


c     cte = -cte_dieter but general dependency in musq.
c   substracted the c_virt
c     -     -7d0 + 1d0/3d0*pi2

      cte=-lns*lns+3d0*lns -7d0! + pi*pi/3d0

c ATTENTION: this last piece pi*pi/3d0 is removed due
c to the change of convention from Gamma(1+eps)-> 1/Gamma[1-eps]

c  First check whether we are diving two zeros. IF so, make substraction
      IF(ABS(RESULTB(1)).LT.1D-4) THEN
      TESTB(1)=ABS(RESULT(1)-CTE*RESULTB(1))
      ELSE
c If not ZERO, then normalized quantities.
      testb(1)=abs(result(1)/(cte*resultb(1))-1d0)
      ENDIF
      
c  First check whether we are diving two zeros. IF so, make substraction
      if(ABS(RESULTB(2)).LT.1D-4) THEN
      TESTB(2)=ABS(RESULT(2)-CTE*RESULTB(2))
      ELSE
c If not ZERO, then normalized quantities.
      testb(2)=abs(result(2)/(cte*resultb(2))-1d0)
      ENDIF
      
      if ( (testb(1).lt.gaugethres) .and. (testb(2).lt.gaugethres)) then
      ward = .true.
      GaugeAccuracy(1)=max(abs(testb(1)),abs(testb(2)))
cFC      print*,'testb(1)', testb(1)
cFC      print*,'testb(2)', testb(2)
c      print*, 'resultb(1)',resultb(1)
c      print*, 'resultb(2)',resultb(2)  
      else
c      print*,'MAXtestb(1)', MAX(testb(1),testb(2))
cFC      print*,'testb(1)', testb(1)
cFC      print*,'testb(2)', testb(2)
cFC      print*, "cte",cte
cFC      print*, 'resultb(1)',resultb(1)
cFC      print*, 'resultb(2)',resultb(2)  
cFC      print*, 'result(1)',result(1)      
cFC      print*, 'result(2)',result(2)      
  
        endif
cFC
cFC      Print*, 'Ward', Ward

      return
      end


      subroutine HexAux(Hex,temp)
      implicit none
      double complex Hex(3)
      double precision temp
      logical ward
      
cFC      print*,"Hex(1)", Hex(1)
cFC      print*,"Hex(2)", Hex(2)
cFC      print*,"Hex(3)", Hex(3)

       ward=.false.

      if(abs(Hex(1)).gt.abs(Hex(2))) then
         if(abs(Hex(1)).gt.abs(Hex(3))) then
            temp=abs((Hex(2)+Hex(3))/(-Hex(1))-1d0)
         else
            temp=abs((Hex(1)+Hex(2))/(-Hex(3))-1d0)
         endif
      else
         if(abs(Hex(2)).gt.abs(Hex(3))) then
            temp=abs((Hex(1)+Hex(3))/(-Hex(2))-1d0)            
         else
            temp=abs((Hex(1)+Hex(2))/(-Hex(3))-1d0)
         endif
      endif

      if(abs(temp).lt.1d-1)then
         ward=.true.
      endif

      end
