      program quadrature1D
      implicit none
      integer::n1,n2
      real::boxLength,mass,valueIntegral
      real,external::pibVal1D
!
      boxLength = float(1)
      mass = float(1)
      write(*,*)' Hrant - x=0.25 ==> ',pibVal1D(1,boxLength,mass,0.25)
      write(*,*)' Hrant - x=0.50 ==> ',pibVal1D(1,boxLength,mass,0.5)
      write(*,*)' Hrant - x=0.75 ==> ',pibVal1D(1,boxLength,mass,0.75)
!
!     Use the trapezoid integration code to numerically integrate Psi*Psi.
!
      call trapezoidInt1D(1,1,boxLength,valueIntegral)
      write(*,*)' Integral(1,1) = ',valueIntegral
      call trapezoidInt1D(1,2,boxLength,valueIntegral)
      write(*,*)' Integral(1,2) = ',valueIntegral
      call trapezoidInt1D(1,12,boxLength,valueIntegral)
      write(*,*)' Integral(1,12) = ',valueIntegral
      call trapezoidInt1D(3,3,boxLength,valueIntegral)
      write(*,*)' Integral(3,3) = ',valueIntegral
      call trapezoidInt1D(3,4,boxLength,valueIntegral)
      write(*,*)' Integral(3,4) = ',valueIntegral
      call trapezoidInt1D(20,20,boxLength,valueIntegral)
      write(*,*)' Integral(20,20) = ',valueIntegral
      call trapezoidInt1D(100,101,boxLength,valueIntegral)
      write(*,*)' Integral(100,101) = ',valueIntegral
!
      end program quadrature1D


      subroutine trapezoidInt1D(n1,n2,l,valueIntegral)
!
      implicit none
      integer,intent(in)::n1,n2
      real,intent(in)::l
      integer::i,nGrid=5
      real::deltaX,xCurrent,valueIntegral
      real,external::pibVal1D
!
 1000 format(3x,'Grid Pt:(,'I3,')  Coord:(',f4.2,')')
      valueIntegral = float(0)
      deltaX = l/nGrid
      xCurrent = float(0)
      do i = 1,nGrid
        write(*,1000) i,xCurrent
        valueIntegral = valueIntegral +  &
          pibVal1D(n1,l,xCurrent)*pibVal1D(n2,l,xCurrent)
        xCurrent = xCurrent + deltaX
        write(*,1000) i,xCurrent
        valueIntegral = valueIntegral +  &
          pibVal1D(n1,l,xCurrent)*pibVal1D(n2,l,xCurrent)
      endDo
      valueIntegral = (deltaX/float(2))*valueIntegral
!
      end subroutine trapezoidInt1D


      function pibVal1D(n,l,x) result(psiX)
!
      implicit none
      integer,intent(in)::n
      real,intent(in)::l,x
      real,intent(out)::psiX
      real::prefactor,temp
      real,save::pi
!
      pi=float(4)*ATan(float(1))
      prefactor = float(2)/l
      prefactor = SQRT(prefactor)
      temp = float(n)*pi*x/l
      psiX = prefactor*Sin(temp)
!
      return
      end
