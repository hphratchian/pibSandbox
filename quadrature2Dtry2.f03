include "mod_quadrature2Dtry2.f03"

      program quadrature2Dtry2
      USE iso_fortran_env
      USE mod_quadrature2Dtry2
      implicit none
      integer(kind=int64)::n1,n2,nGrid
      integer(kind=int64)::IntegrationSegmentsX,IntegrationSegmentsY
      real(kind=real64)::boxLength,mass,valueIntegral
      real(kind=real64)::x
      character(len=256)::commandLineArg
      logical::Debug=.false.
!
 2000 format(1x,'Integral(',I4,',',I4,') = ',F12.5)
 8000 format(3x,I5,',',I5,':   ',F15.8)
!
!     Hardwire the box length and mass.
!
      boxLength = float(1)
      mass = float(1)
!
!     If requested by the DEBUG flag, print out some values for the 1D PID
!     wave function code.
!
      if(DEBUG) then
        x = 0.25
        write(*,*)' Hrant - x=0.25 ==> ',pibVal1D(1_int64,boxLength,x)
        x = 0.5
        write(*,*)' Hrant - x=0.50 ==> ',pibVal1D(1_int64,boxLength,x)
        x = 0.75
        write(*,*)' Hrant - x=0.75 ==> ',pibVal1D(1_int64,boxLength,x)
      endIf
!
!     If requested by the DEBUG flag, run the original trapezoid integration
!     code to numerically integrate Psi*Psi.
!
      if(DEBUG) then
        call trapezoidInt2D(1_int64,1_int64,1_int64,boxLength,mass,valueIntegral)
        write(*,2000) 1,1,valueIntegral
        call trapezoidInt2D(1_int64,1_int64,2_int64,boxLength,mass,valueIntegral)
        write(*,2000) 1,2,valueIntegral
        call trapezoidInt2D(1_int64,2_int64,2_int64,boxLength,mass,valueIntegral)
        write(*,2000) 2,2,valueIntegral
        call trapezoidInt2D(1_int64,1_int64,3_int64,boxLength,mass,valueIntegral)
        write(*,2000) 1,3,valueIntegral
      endIf
      call trapezoidInt2D(1_int64,2_int64,2_int64,boxLength,mass,valueIntegral)
      write(*,2000) 2,2,valueIntegral
!
!     Try out the new quadrature code...
!
      call GET_COMMAND_ARGUMENT(1,commandLineArg)
      read(commandLineArg,*) IntegrationSegmentsX
      call GET_COMMAND_ARGUMENT(2,commandLineArg)
      read(commandLineArg,*) IntegrationSegmentsY
      call quadratureInitialize(2_int64)
      call quadratureSetGrid_trapezoidInt2D(IntegrationSegmentsX,  &
        IntegrationSegmentsY,0.0_real64,0.0_real64,boxLength,boxLength,  &
        DebugPrint=.false.)
      call quadratureGeneral(valueIntegral,DebugPrint=.false.)
      write(*,8000) IntegrationSegmentsX,IntegrationSegmentsY,valueIntegral
!
      end program quadrature2Dtry2
