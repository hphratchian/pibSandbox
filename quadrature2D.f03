      program quadrature2D
      USE iso_fortran_env
      implicit none
      integer::n1,n2
      real(kind=real64)::boxLength,mass,valueIntegral
      real(kind=real64)::x
      real(kind=real64),external::pibVal1D
!
 2000 format(1x,'Integral(',I4,',',I4,') = ',F12.5)
!
      boxLength = float(1)
      mass = float(1)
      x = 0.25
      write(*,*)' Hrant - x=0.25 ==> ',pibVal1D(1,boxLength,x)
      x = 0.5
      write(*,*)' Hrant - x=0.50 ==> ',pibVal1D(1,boxLength,x)
      x = 0.75
      write(*,*)' Hrant - x=0.75 ==> ',pibVal1D(1,boxLength,x)
!
!     Use the trapezoid integration code to numerically integrate Psi*Psi.
!
      call trapezoidInt2D(1,1,1,boxLength,mass,valueIntegral)
      write(*,2000) 1,1,valueIntegral
      call trapezoidInt2D(1,1,2,boxLength,mass,valueIntegral)
      write(*,2000) 1,2,valueIntegral
      call trapezoidInt2D(1,2,2,boxLength,mass,valueIntegral)
      write(*,2000) 2,2,valueIntegral
      call trapezoidInt2D(1,1,3,boxLength,mass,valueIntegral)
      write(*,2000) 1,3,valueIntegral

!hph+
!      call trapezoidInt2D(1,1,12,boxLength,mass,valueIntegral)
!      write(*,*)' Integral(1,12) = ',valueIntegral
!      call trapezoidInt2D(1,3,3,boxLength,mass,valueIntegral)
!      write(*,*)' Integral(3,3) = ',valueIntegral
!      call trapezoidInt2D(1,3,4,boxLength,mass,valueIntegral)
!      write(*,*)' Integral(3,4) = ',valueIntegral
!      call trapezoidInt2D(1,20,20,boxLength,mass,valueIntegral)
!      write(*,*)' Integral(20,20) = ',valueIntegral
!      call trapezoidInt2D(1,100,101,boxLength,mass,valueIntegral)
!      write(*,*)' Integral(100,101) = ',valueIntegral
!hph-

!
      end program quadrature2D


      subroutine trapezoidInt2D(IntegralType,n1,n2,l,m,valueIntegral,  &
        nGridInX,nGridInY)
      USE iso_fortran_env
!
!     The argument <IntegralType> is a flag indicated the type of integral that
!     should be integrated here. The available options are:
!           1. Overlap Integral
!           2. Diagonal 2-Center Coulomb Integral
!
!
      implicit none
      integer,intent(in)::IntegralType,n1,n2
      integer,intent(in),OPTIONAL::nGridInX,nGridInY
      real(kind=real64),intent(in)::l,m
      integer::iX,iY
      integer::nGridX,nGridY
      integer,parameter::nGridMultiplier=1000
      real(kind=real64)::deltaX,deltaY,xCurrent,yCurrent,xNext,yNext,  &
        valueR12,valueVolume11,valueVolume21,valueVolume12,valueVolume22,  &
        valueIntegral
      real(kind=real64),parameter::small=float(1)/float(10000)
      real(kind=real64),external::pibVal1D
      real(kind=real64),external::pib2DIntegrand
      logical::DEBUG=.false.
!
 1000 format(/,1x,'Enter trapezoidInt2D')
 1100 format(2x,'nGridX=',I5,3x,'nGridY=',I5,/,3x,'deltaX=',ES10.3,3x,  &
   'deltaY=',ES10.3)
 2000 format(3x,'Grid Pt:(,'I3,',',I3,')  11 Coord:(',f4.2,',',f4.2,')')
 2010 format(23x,I2,' Coord:(',F4.2,',',F4.2,')',)
 8000 format(1x,'*** WARNING ***',/,3x,'Routine trapezoidInt2D')
!
      if(DEBUG) write(*,1000)
      valueIntegral = float(0)
      if(Present(nGridInX)) then
        nGridX = nGridInX
      else
        nGridX = nGridMultiplier*2*n1
      endIf
      if(Present(nGridInY)) then
        nGridY = nGridInY
      else
        nGridY = nGridMultiplier*2*n2
      endIf
      deltaX = l/float(nGridX)
      deltaY = l/float(nGridY)
      if(DEBUG) write(*,1100) nGridX,nGridY,deltaX,deltaY
      if(deltaX.le.float(10)*small) then
        write(*,8000)
        write(*,1100) nGridX,nGridY,deltaX,deltaY
      endIf
      xCurrent = float(0)
      yCurrent = float(0)
      xNext = xCurrent + deltaX
      yNext = yCurrent + deltaY
      do iX = 1,nGridX
        do iY = 1,nGridY
          if(DEBUG) then
            write(*,2000) nGridX,nGridY,xCurrent,yCurrent
            write(*,2010) 21,xNext,yCurrent
            write(*,2010) 21,xCurrent,yNext
            write(*,2010) 22,xNext,yNext
          endIf
          valueVolume11 = pib2DIntegrand(n1,n2,l,l,xCurrent,yCurrent)
          valueVolume21 = pib2DIntegrand(n1,n2,l,l,xNext,yCurrent)
          valueVolume12 = pib2DIntegrand(n1,n2,l,l,xCurrent,yNext)
          valueVolume22 = pib2DIntegrand(n1,n2,l,l,xNext,yNext)
          valueIntegral = valueIntegral + valueVolume11 + valueVolume21  &
            + valueVolume12 + valueVolume22
          yCurrent = yNext
          yNext = yCurrent + deltaY
        endDo
        xCurrent = xNext
        yCurrent = float(0)
        xNext = xCurrent + deltaX
        yNext = yCurrent + deltaY
      endDo
      valueIntegral = (deltaX*deltaY/float(4))*valueIntegral
!
      end subroutine trapezoidInt2D


      function pib2DIntegrand(nx,ny,lx,ly,x,y) result(valueIntegrand)
      USE iso_fortran_env
!
      implicit none
      integer,intent(in)::nx,ny
      real(kind=real64),intent(in)::lx,ly,x,y
      real(kind=real64),intent(out)::valueIntegrand
      real(kind=real64),external::pibVal1D
!
      valueIntegrand = pibVal1D(nx,lx,x)*pibVal1D(ny,ly,y)
!
      return
      end function pib2DIntegrand


      function pibVal2D(nx,ny,lx,ly,m,x,y) result(psiXY)
      USE iso_fortran_env
!
      implicit none
      integer,intent(in)::nx,ny
      real(kind=real64),intent(in)::lx,ly,m,x,y
      real(kind=real64),intent(out)::psiXY
      real(kind=real64),external::pibVal1D
!
      psiXY = pibVal1D(nx,lx,x)*pibVal1D(ny,ly,y)
!
      return
      end function pibVal2D


      function pibVal1D(n,l,x) result(psiX)
      USE iso_fortran_env
!
      implicit none
      integer,intent(in)::n
      real(kind=real64),intent(in)::l,x
      real(kind=real64),intent(out)::psiX
      real(kind=real64)::prefactor,temp
      real(kind=real64),save::pi
!
      pi=float(4)*ATan(float(1))
      prefactor = float(2)/l
      prefactor = SQRT(prefactor)
      temp = float(n)*pi*x/l
      psiX = prefactor*Sin(temp)
!
      return
      end function pibVal1D
