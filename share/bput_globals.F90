  ! this module contains necessary methods and data sturctures that
  ! implement pnet_cdf bput.
  !
  ! all methods's names begin with "Bput"
  module bput_globals
      implicit none

      integer, allocatable :: reqs_(:)
      integer :: num_reqs_
      integer :: i_
      integer :: allocated_
      contains

      subroutine BputInitGlobals(num_reqs)
          ! init reqs_ to have size num_reqs; possibly all -1?
          integer, intent(in) :: num_reqs
          
          if (allocated(reqs_)) deallocate(reqs_)
          allocate(reqs_(num_reqs))

          num_reqs_ = num_reqs
          i_ = 1
      end subroutine BputInitGlobals

      subroutine BputSetNextReqVal(req_val)
          integer, intent(in) :: req_val

          if (i_ .LE. num_reqs_) then
              reqs_(i_) = req_val
              i_ = i_ + 1
          ! else
              ! error?
          endif
          
      end subroutine BputSetNextReqVal

      subroutine BputGetNumOfReqs(num_reqs)
          integer, intent(out) :: num_reqs

          num_reqs = num_reqs_
      end subroutine BputGetNumOfReqs

      subroutine BputGetReqVal(req_idx, req_val)
          integer, intent(in) :: req_idx
          integer, intent(out) :: req_val

          if (req_idx .LE. num_reqs_) then
              req_val = reqs_(req_idx)
          else
              ! error?
              req_val = -1
          endif
      end subroutine BputGetReqVal

      subroutine BputLowerCase(MemoryOrder,MemOrd)
          character*(*) ,intent(in)  :: MemoryOrder
          character*(*) ,intent(out) :: MemOrd
          character*1                :: c
          integer       ,parameter   :: upper_to_lower =IACHAR('a')-IACHAR('A')
          integer                    :: i,N
          
          MemOrd = ' '
          N = len(MemoryOrder)
          MemOrd(1:N) = MemoryOrder(1:N)
          do i=1,N
              c = MemoryOrder(i:i)
              if('A'<=c .and. c <='Z') MemOrd(i:i)=achar(iachar(c)+upper_to_lower)
          enddo
          return
      end subroutine BputLowerCase

      SUBROUTINE BputGetDim(MemoryOrder,NDim,Status)
          include 'wrf_status_codes.h'
          character*(*) ,intent(in)  :: MemoryOrder
          integer       ,intent(out) :: NDim
          integer       ,intent(out) :: Status
          character*3                :: MemOrd
        
          call BputLowerCase(MemoryOrder,MemOrd)
          select case (MemOrd)
            case ('xyz','xzy','yxz','yzx','zxy','zyx','xsz','xez','ysz','yez')
              NDim = 3
            case ('xy','yx','xs','xe','ys','ye')
              NDim = 2
            case ('z','c')
              NDim = 1
            case ('0')  ! NDim=0 for scalars.  TBH:  20060502
              NDim = 0
            case default
              print *, 'memory order = ',MemOrd,'  ',MemoryOrder
              Status = WRF_WARN_BAD_MEMORYORDER
              return
          end select
          Status = WRF_NO_ERR
          return
      END SUBROUTINE BputGetDim

      SUBROUTINE BputGetGridSize(p, sizeOut)
          USE module_domain_type, ONLY : fieldlist
          TYPE(fieldlist), intent(inout), POINTER :: p
          integer, intent(out) :: sizeOut
          INTEGER :: ndim, ierr

          CALL BputGetDim(TRIM(p%MemoryOrder), ndim, ierr)
      
          sizeOut = 0
          IF (ndim .EQ. 0) THEN    
            sizeOut = 1
          ELSE IF (ndim .EQ. 1) THEN
            if (p%ep1 - p%sp1 > 0) &
            sizeOut = (p%ep1 - p%sp1 + 1)
          ELSE IF (ndim .EQ. 2) THEN
            if (p%ep1 - p%sp1 > 0 .AND. p%ep2 - p%sp2 > 0) &
            sizeOut = (p%ep1 - p%sp1 + 1) * (p%ep2 - p%sp2 + 1) 
          ELSE IF (ndim .EQ. 3) THEN
            if (p%ep1 - p%sp1 > 0 .AND. p%ep2 - p%sp2 > 0 .AND. p%ep3 - p%sp3 > 0) &
            sizeOut = (p%ep1 - p%sp1 + 1) * (p%ep2 - p%sp2 + 1) * (p%ep3 - p%sp3 + 1)
          ENDIF
          return
      END SUBROUTINE BputGetGridSize

      SUBROUTINE BputGetSizeOfType(Type, Size)
          CHARACTER*1, intent(in) :: Type
          INTEGER, intent(out) :: Size
          ! LOGICAL :: temp_logical
          REAL :: temp_real
          REAL*8 :: temp_double
          
          Size = 1
          IF (Type .EQ. 'r') THEN
            Size = SIZEOF(temp_real)
          !   CALL test_called(Size)
          ELSE IF (Type .EQ. 'd') THEN
            Size = SIZEOF(temp_double)
          ELSE IF (Type .EQ. 'i') THEN
            Size = SIZEOF(Size)
          ELSE IF (Type .EQ. 'l') THEN
            Size = SIZEOF(Size)
          ENDIF
          RETURN
      END SUBROUTINE BputGetSizeOfType

      subroutine BputGetBufferSizeAndNumCalls(fieldListPtr,switch,grid, totalSizeOut, numCallsOut)
          USE module_domain_type, ONLY : fieldlist
          USE module_domain
          USE module_io

          ! input
          TYPE(fieldlist), intent(in), POINTER :: fieldListPtr
          INTEGER, INTENT(IN) :: switch
          TYPE(domain), intent(in) :: grid

          ! output
          integer, intent(out):: totalSizeOut
          integer, intent(out):: numCallsOut

          integer :: gridSize, typeSize
          INTEGER newSwitch, itrace
          TYPE(fieldlist), POINTER :: p
          p => fieldListPtr
          
          ! init variables
          newSwitch = switch
          totalSizeOut = 0
          numCallsOut = 0

          DO WHILE ( ASSOCIATED( p ) )
            IF ( p%ProcOrient .NE. 'X' .AND. p%ProcOrient .NE. 'Y' ) THEN
              gridSize = 0
              IF (p%Ndim .EQ. 0)  THEN
                ! TODO: check 0 or 1
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    ! set grid szie. equals 1 when p%Ndim == 0
                    gridSize = 1
                    
                    typeSize = 4
                    CALL BputGetSizeOfType(p%Type, typeSize)
                    
                    totalSizeOut = totalSizeOut + gridSize * typeSize
                    numCallsOut = numCallsOut + 1
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 1) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      CALL BputGetGridSize(p, gridSize)

                      typeSize = 4
                      CALL BputGetSizeOfType(p%Type, typeSize)
                  
                      totalSizeOut = totalSizeOut + gridSize * typeSize
                      numCallsOut = numCallsOut + 1
                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 2) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) .AND.  &
                    ( .NOT. p%subgrid_x .OR. (p%subgrid_x .AND. grid%sr_x .GT. 0) ) .AND. &
                    ( .NOT. p%subgrid_y .OR. (p%subgrid_y .AND. grid%sr_y .GT. 0) )       &
                  ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      CALL BputGetGridSize(p, gridSize)

                      typeSize = 4
                      CALL BputGetSizeOfType(p%Type, typeSize)
                      totalSizeOut = totalSizeOut + gridSize * typeSize
                      numCallsOut = numCallsOut + 1

                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 3) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) .AND.  &
                    ( .NOT. p%subgrid_x .OR. (p%subgrid_x .AND. grid%sr_x .GT. 0) ) .AND. &
                    ( .NOT. p%subgrid_y .OR. (p%subgrid_y .AND. grid%sr_y .GT. 0) )       &
                  ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      CALL BputGetGridSize(p, gridSize)

                      typeSize = 4
                      CALL BputGetSizeOfType(p%Type, typeSize)
                      totalSizeOut = totalSizeOut + gridSize * typeSize
                      numCallsOut = numCallsOut + 1
                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 4 .AND. p%scalar_array ) THEN
                IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                  DO itrace = PARAM_FIRST_SCALAR , p%num_table(grid%id)
                    
                    IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams_table(grid%id,itrace)%stream,newSwitch)) THEN
                      CALL BputGetGridSize(p, gridSize)

                      typeSize = 4
                      CALL BputGetSizeOfType(p%Type, typeSize)
                      totalSizeOut = totalSizeOut + gridSize * typeSize
                      numCallsOut = numCallsOut + 1
                    ENDIF
                  ENDDO
                  
                ENDIF
              ENDIF

            ENDIF
            p => p%next
          ENDDO
        
      end subroutine BputGetBufferSizeAndNumCalls

  end module