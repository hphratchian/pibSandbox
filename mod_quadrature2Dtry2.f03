      module mod_quadrature2Dtry2
      USE iso_fortran_env
      implicit none
      real(kind=real64)::quadtratureDimension,quadratureNGridPoints
      real(kind=real64),dimension(:),allocatable::quadratureWeights,quadratureGridSize
      real(kind=real64),dimension(:,:),allocatable::quadratureGrid
      integer(kind=int64),dimension(:,:),allocatable::quadratureGridLabels

      
      CONTAINS
!
!
!
!PROCEDURE quadratureInitialize
      subroutine quadratureInitialize(nDimensions)
!
!     This subroutine is used to initialize the quadrature object before
!     determining the grid points and the weights.
!
      implicit none
      integer(kind=int64),intent(in)::nDimensions
!
!     Begin by setting the object's dimension flag equal to the value in
!     nDimensions sent to the routine. Then, dealocate the arrays in the object
!     if they've been previously allocated.
!
      quadtratureDimension = nDimensions
      if(Allocated(quadratureWeights)) DeAllocate(quadratureWeights)
      if(Allocated(quadratureGridSize)) DeAllocate(quadratureGridSize)
      if(Allocated(quadratureGrid)) DeAllocate(quadratureGrid)
      if(Allocated(quadratureGridLabels)) DeAllocate(quadratureGridLabels)
!
!     Allocate memory for quadratureGridSize in the object.
!
      Allocate(quadratureGridSize(nDimensions))
!
      return
      end subroutine quadratureInitialize

!
!PROCEDURE quadratureSetGrid_trapezoidInt2D
      subroutine quadratureSetGrid_trapezoidInt2D(nIntegrationSegmentsX,  &
        nIntegrationSegmentsY,x0,y0,x1,y1,DebugPrint)
!
!     This subroutine sets up the quadrature grid points and integration weights
!     for the 2D trapezoid method. As input, the calling program unit sends the
!     number of integration segments in each dimension (labeled as X and Y) and
!     the limits of integration in each dimension.
!
      implicit none
      integer(kind=int64),intent(in)::nIntegrationSegmentsX,  &
        nIntegrationSegmentsY
      real(kind=real64),intent(in)::x0,y0,x1,y1
      logical,optional,intent(in)::DebugPrint
      integer(kind=int64)::nGridPointsX,nGridPointsY,iX,iY,iTypeX,  &
        iTypeY,iGridPoint
      real(kind=real64)::stepSizeX,stepSizeY,xVal,yVal,area
      logical::debug=.false.
!
 2000 format(3x,'Grid Point No. ',I6,'  (',I6,',',I6,') ',&
        '==> (',F6.3,',',F6.3,')  |  w = ',F6.3,'.')
!
!     Set-up debug printing flag.
!
      if(PRESENT(DebugPrint)) DEBUG = DebugPrint
!
!     Begin by setting the elements of quadratureGridSize in the quadrature
!     object and then allocating arrays for the grid and weights in the
!     quadrature object.
!
      nGridPointsX = nIntegrationSegmentsX + 1
      nGridPointsY = nIntegrationSegmentsY + 1
      quadratureGridSize(1) = nGridPointsX
      quadratureGridSize(2) = nGridPointsY
      quadratureNGridPoints = nGridPointsX*nGridPointsY
      Allocate(quadratureWeights(nGridPointsX*nGridPointsY))
      Allocate(quadratureGrid(nGridPointsX*nGridPointsY,2),  &
        quadratureGridLabels(nGridPointsX*nGridPointsY,2))
!
!     Determine the integration step size in each dimension and evaluate the
!     area value that will be used in evaluating the integration weights.
!
      stepSizeX = (x1-x0)/float(nIntegrationSegmentsX)
      stepSizeY = (y1-y0)/float(nIntegrationSegmentsY)
      area = stepSizeX*stepSizeY
!
!     Loop through the grid points and fill out the quadrature info into the
!     object.
!
      iGridPoint = 0
      do iX = 0,nGridPointsX-1
        do iY = 0,nGridPointsY-1
          iGridPoint = iGridPoint + 1
          xVal = x0 + stepSizeX*iX
          yVal = y0 + stepSizeY*iY
          quadratureGrid(iGridPoint,1) = xVal
          quadratureGrid(iGridPoint,2) = yVal
          iTypeX = 0
          if(iX.ne.0.and.iX.ne.nGridPointsX-1) iTypeX = 1
          iTypeY = 0
          if(iY.ne.0.and.iY.ne.nGridPointsY-1) iTypeY = 1
          select case(iTypeX)
          case (0)
            select case(iTypeY)
            case(0)
              quadratureWeights(iGridPoint) = area*dfloat(1)/dfloat(4)
            case default
              quadratureWeights(iGridPoint) = area*dfloat(1)/dfloat(2)
            end select
          case default
            select case(iTypeY)
            case(0)
              quadratureWeights(iGridPoint) = area*dfloat(1)/dfloat(2)
            case default
              quadratureWeights(iGridPoint) = area*dfloat(1)
            end select
          end select
          if(DEBUG) write(*,2000) iGridPoint,iX,iY,xVal,yVal,  &
            quadratureWeights(iGridPoint)
        endDo
      endDo
!
      end subroutine quadratureSetGrid_trapezoidInt2D

!
!PROCEDURE quadratureGeneral
!hph      function quadratureGeneral(integrandFunction) result(quadratureVal)
!hph      function quadratureGeneral() result(quadratureVal)
      subroutine quadratureGeneral(quadratureVal,DebugPrint)
!
!     This function evaluates an integral using the general quadrature objects.
!     The function dummy argument integrandFunction is sent for the integrand.
!
!
      implicit none
!hph      real(kind=real64),external::integrandFunction
      real(kind=real64)::quadratureVal
      logical,optional,intent(in)::DebugPrint
      integer(kind=int64)::i
      real(kind=real64),dimension(:),allocatable::gridValues
      logical::DEBUG=.false.
!
 1000 format(3x,i6,': (',f6.3,',',f6.3,') = ',f12.8)
!
!     Set-up debug printing flag.
!
      if(PRESENT(DebugPrint)) DEBUG = DebugPrint
!
!     Build the integrand values at the grid points. Then, take the dot product
!     of with the quadrature weights to get the final integral value that's put
!     into quadratureVal.
!
      Allocate(gridValues(quadratureNGridPoints))
      do i = 1,quadratureNGridPoints
        gridValues(i) = pib2DIntegrand(2_int64,2_int64,1.0_real64,  &
          1.0_real64,quadratureGrid(i,1),quadratureGrid(i,2))
        if(DEBUG)  &
          write(*,1000) i,quadratureGrid(i,1),quadratureGrid(i,2),gridValues(i)
      endDo
      quadratureVal = dot_product(quadratureWeights,gridValues)
      DeAllocate(gridValues)
!
      return
      end subroutine quadratureGeneral
!hph      end function quadratureGeneral

!
!PROCEDURE trapezoidInt2D
      subroutine trapezoidInt2D(IntegralType,n1,n2,l,m,valueIntegral,  &
        nGridInX,nGridInY,DebugPrint)
!
!     The argument <IntegralType> is a flag indicated the type of integral that
!     should be integrated here. The available options are:
!           1. Overlap Integral
!           2. Diagonal 2-Center Coulomb Integral
!
!
      implicit none
      integer(kind=int64),intent(in)::IntegralType,n1,n2
      integer(kind=int64),intent(in),OPTIONAL::nGridInX,nGridInY
      real(kind=real64),intent(in)::l,m
      logical,optional,intent(in)::DebugPrint
      integer::iX,iY
      integer::nGridX,nGridY
      integer(kind=int64),parameter::nGridMultiplier=1000
      real(kind=real64)::deltaX,deltaY,xCurrent,yCurrent,xNext,yNext,  &
        valueR12,valueVolume11,valueVolume21,valueVolume12,valueVolume22,  &
        valueIntegral
      real(kind=real64),parameter::small=float(1)/float(10000)
      logical::DEBUG=.false.
!
 1000 format(/,1x,'Enter trapezoidInt2D')
 1100 format(2x,'nGridX=',I5,3x,'nGridY=',I5,/,3x,'deltaX=',ES10.3,3x,  &
   'deltaY=',ES10.3)
 2000 format(3x,'Grid Pt:(,'I3,',',I3,')  11 Coord:(',f4.2,',',f4.2,')')
 2010 format(23x,I2,' Coord:(',F4.2,',',F4.2,')',)
 8000 format(1x,'*** WARNING ***',/,3x,'Routine trapezoidInt2D')
!
!     Set-up debug printing flag.
!
      if(PRESENT(DebugPrint)) DEBUG = DebugPrint
!
!     Figure out some of the key terms used below.
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
      write(*,1100) nGridX,nGridY,deltaX,deltaY
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

!
!PROCEDURE pib2DIntegrandWrapper
      function pib2DIntegrandWrapper(coordinates,integerFlags,realFlags)  &
        result(valueIntegrand)
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::coordinates
      integer(kind=int64),dimension(:),intent(inOut)::integerFlags
      real(kind=real64),dimension(:),intent(inOut)::realFlags
      real(kind=real64),intent(out)::valueIntegrand
!
      valueIntegrand = pib2DIntegrand(integerFlags(1),integerFlags(2),  &
        realFlags(1),realFlags(2),coordinates(1),coordinates(2))
!
      return
      end function pib2DIntegrandWrapper

!
!PROCEDURE pib2DIntegrand
      function pib2DIntegrand(nx,ny,lx,ly,x,y) result(valueIntegrand)
!
      implicit none
      integer(kind=int64),intent(in)::nx,ny
      real(kind=real64),intent(in)::lx,ly,x,y
      real(kind=real64),intent(out)::valueIntegrand
!
      valueIntegrand = pibVal1D(nx,lx,x)*pibVal1D(ny,ly,y)
!
      return
      end function pib2DIntegrand

!
!PROCEDURE pibVal2D
      function pibVal2D(nx,ny,lx,ly,m,x,y) result(psiXY)
!
      implicit none
      integer(kind=int64),intent(in)::nx,ny
      real(kind=real64),intent(in)::lx,ly,m,x,y
      real(kind=real64),intent(out)::psiXY
!
      psiXY = pibVal1D(nx,lx,x)*pibVal1D(ny,ly,y)
!
      return
      end function pibVal2D

!
!PROCEDURE pibVal1D
      function pibVal1D(n,l,x) result(psiX)
!
      implicit none
      integer(kind=int64),intent(in)::n
      real(kind=real64),intent(in)::l,x
      real(kind=real64),intent(out)::psiX
      real(kind=real64)::prefactor,temp
      real(kind=real64),save::pi
!
      if(n.le.0) then
        psiX = 1.0_real64
        return
      endIf
      pi=float(4)*ATan(float(1))
      prefactor = float(2)/l
      prefactor = SQRT(prefactor)
      temp = float(n)*pi*x/l
      psiX = prefactor*Sin(temp)
!
      return
      end function pibVal1D

      end module mod_quadrature2Dtry2
