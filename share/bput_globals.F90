  ! this module contains necessary methods and data sturctures that
  ! implement pnet_cdf non-blocking bput.
  !
  ! all methods's names begin with "Bput"

  module bput_globals
      implicit none

      integer, allocatable :: reqs_(:)
      integer :: num_reqs_
      integer :: i_
      integer :: use_bput_  ! 1 use bput; 0 do not use bput (e.g. dryrun)
      CHARACTER*256 err_msg
      integer, parameter :: DateStrLen = 19
      contains

      subroutine BputAttach(fid, bput_buffer_size, err)
        use pnetcdf
        integer, intent(in) :: fid
        integer*8, intent(in) :: bput_buffer_size
        integer, intent(out):: err
        err = nfmpi_buffer_attach(fid, bput_buffer_size)
      end subroutine BputAttach

      subroutine BputDetach(fid)
        use pnetcdf
        integer, intent(in) :: fid
        integer:: err
        err = nfmpi_buffer_detach(fid)
      end subroutine BputDetach

      subroutine BputWaitAll(ncid)
        use pnetcdf
        integer, intent(in) :: ncid
        integer, allocatable :: st(:)
        integer :: err
        allocate(st(num_reqs_))
        err = nf90mpi_wait_all(ncid, num_reqs_, reqs_, st)
        deallocate(st)
      end subroutine BputWaitAll


      subroutine BputSetUse(use_bput)
        integer, intent(in) :: use_bput
        use_bput_ = use_bput
      end subroutine BputSetUse

      subroutine BputGetUse(use_bput_out)
        integer, intent(out) :: use_bput_out
        use_bput_out = use_bput_
      end subroutine BputGetUse

      subroutine BputSetNumReqs(num_reqs)
          integer, intent(in) :: num_reqs
          
          if (allocated(reqs_)) deallocate(reqs_)
          allocate(reqs_(num_reqs))

          num_reqs_ = num_reqs
          i_ = 1
      end subroutine BputSetNumReqs

      subroutine BputCleanGlobals()
        if (allocated(reqs_)) deallocate(reqs_)
        num_reqs_ = 0
        i_ = 1
        use_bput_ = 0
      end subroutine BputCleanGlobals


      subroutine BputSetNextReqVal(req_val)
          integer, intent(in) :: req_val

          if (i_ .LE. num_reqs_) then
              reqs_(i_) = req_val
              i_ = i_ + 1
          endif
          
      end subroutine BputSetNextReqVal

      subroutine BputGetNumOfReqs(num_reqs)
          integer, intent(out) :: num_reqs

          num_reqs = num_reqs_
      end subroutine BputGetNumOfReqs

      subroutine BputGetAllReqs(all_reqs)
        integer, intent(inout), allocatable :: all_reqs(:)
        integer :: i
        if (allocated(all_reqs)) deallocate(all_reqs)
        allocate(all_reqs(num_reqs_))

        do i = 1 , num_reqs_
          all_reqs(i) = reqs_(i)
        enddo

      end subroutine BputGetAllReqs

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
              Status = -1
              return
          end select
          Status = 0
          return
      END SUBROUTINE BputGetDim

      SUBROUTINE BputGetGridSize(p, sizeOut)
          USE module_domain_type, ONLY : fieldlist
          TYPE(fieldlist), intent(inout), POINTER :: p
          integer, intent(out) :: sizeOut
          INTEGER :: ndim, ierr

          CALL BputGetDim(TRIM(p%MemoryOrder), ndim, ierr)
          IF (ierr < 0) THEN
            sizeOut = -1
            return
          ENDIF
      
          sizeOut = 0
          IF (ndim .EQ. 0) THEN    
            sizeOut = 1
          ELSE IF (ndim .EQ. 1) THEN
            if (p%ep1 - p%sp1 >= 0) &
            sizeOut = (p%ep1 - p%sp1 + 1)
          ELSE IF (ndim .EQ. 2) THEN
            if (p%ep1 - p%sp1 >= 0 .AND. p%ep2 - p%sp2 >= 0) &
            sizeOut = (p%ep1 - p%sp1 + 1) * (p%ep2 - p%sp2 + 1) 
          ELSE IF (ndim .EQ. 3) THEN
            if (p%ep1 - p%sp1 >= 0 .AND. p%ep2 - p%sp2 >= 0 .AND. p%ep3 - p%sp3 >= 0) &
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

          ! input variables
          TYPE(fieldlist), intent(in), POINTER :: fieldListPtr
          INTEGER, INTENT(IN) :: switch
          TYPE(domain), intent(in) :: grid

          ! output variables
          integer*8, intent(out):: totalSizeOut
          integer, intent(out):: numCallsOut

          ! temp variables
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
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    gridSize = 1 ! grid szie=1 when p%Ndim == 0
                    typeSize = 4 ! init default type size
                    CALL BputGetSizeOfType(p%Type, typeSize)

                    totalSizeOut = totalSizeOut + gridSize * typeSize
                    numCallsOut = numCallsOut + 1
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 1) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newSwitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      CALL BputGetGridSize(p, gridSize) ! get grid size
                      IF (gridSize >= 0) THEN
                        typeSize = 4 ! init default type size
                        CALL BputGetSizeOfType(p%Type, typeSize)

                        totalSizeOut = totalSizeOut + gridSize * typeSize
                        numCallsOut = numCallsOut + 1
                      ENDIF
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
                      CALL BputGetGridSize(p, gridSize) ! get grid size
                      IF (gridSize >= 0) THEN
                        typeSize = 4 ! init default type size
                        CALL BputGetSizeOfType(p%Type, typeSize)

                        totalSizeOut = totalSizeOut + gridSize * typeSize
                        numCallsOut = numCallsOut + 1
                      ENDIF
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
                      CALL BputGetGridSize(p, gridSize) ! get grid size
                      IF (gridSize >= 0) THEN
                        typeSize = 4 ! init default type size
                        CALL BputGetSizeOfType(p%Type, typeSize)

                        totalSizeOut = totalSizeOut + gridSize * typeSize
                        numCallsOut = numCallsOut + 1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 4 .AND. p%scalar_array ) THEN
                IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                  DO itrace = PARAM_FIRST_SCALAR , p%num_table(grid%id)
                    
                    IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams_table(grid%id,itrace)%stream,newSwitch)) THEN
                      CALL BputGetGridSize(p, gridSize) ! get grid size
                      IF (gridSize >= 0) THEN
                        typeSize = 4 ! init default type size
                        CALL BputGetSizeOfType(p%Type, typeSize)

                        totalSizeOut = totalSizeOut + gridSize * typeSize
                        numCallsOut = numCallsOut + 1
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
            p => p%next
          ENDDO

          ! DateStr
          totalSizeOut = totalSizeOut + DateStrLen
          numCallsOut = numCallsOut + 1
        
      end subroutine BputGetBufferSizeAndNumCalls

      subroutine BputGetNCID(DataHandle,NcidOut)
        use wrf_data_pnc
        USE module_io
        include 'wrf_status_codes.h'
        integer               ,intent(in)     :: DataHandle
        type(wrf_data_handle) ,pointer        :: DH
        integer, intent(out) :: NcidOut
        integer :: hndl

        if ( DataHandle .GE. 1 .AND. DataHandle .LE. 1000 ) then
          hndl = wrf_io_handles(DataHandle) ! get true handle (integer)
        else
          hndl = -1
          NcidOut = -1
          return
        endif

        if(hndl < 1 .or. hndl > WrfDataHandleMax) then
          NcidOut = -1
          return
        endif

        DH => WrfDataHandles(hndl)
        if(DH%Free) then
          NcidOut = -1
          return
        endif

        NcidOut = DH%NCID
        return
      end subroutine BputGetNCID
  end module