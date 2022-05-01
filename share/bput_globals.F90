  module bput_globals
      implicit none

      integer, allocatable :: reqs_(:)
      integer :: num_reqs_
      integer :: i_
      integer :: allocated_
      contains

      subroutine init_bput_globals(num_reqs)
          ! init reqs_ to have size num_reqs; possibly all -1?
          integer, intent(in) :: num_reqs
          
          if (allocated(reqs_)) deallocate(reqs_)
          allocate(reqs_(num_reqs))

          num_reqs_ = num_reqs
          i_ = 1
      end subroutine init_bput_globals

      subroutine set_bput_globals_req(req_val)
          integer, intent(in) :: req_val

          if (i_ .LE. num_reqs_) then
              reqs_(i_) = req_val
              i_ = i_ + 1
          ! else
              ! error?
          endif
          
      end subroutine set_bput_globals_req

      subroutine get_bput_globals_num_reqs(num_reqs)
          integer, intent(out) :: num_reqs

          num_reqs = num_reqs_
      end subroutine get_bput_globals_num_reqs

      subroutine get_bput_globals_reqs_at_i(req_idx, req_val)
          integer, intent(in) :: req_idx
          integer, intent(out) :: req_val

          if (req_idx .LE. num_reqs_) then
              req_val = reqs_(req_idx)
          else
              ! error?
              req_val = -1
          endif
      end subroutine get_bput_globals_reqs_at_i

      subroutine ZanhuaLowerCase(MemoryOrder,MemOrd)
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
      end subroutine ZanhuaLowerCase

      SUBROUTINE ZanhuaGetDim(MemoryOrder,NDim,Status)
          include 'wrf_status_codes.h'
          character*(*) ,intent(in)  :: MemoryOrder
          integer       ,intent(out) :: NDim
          integer       ,intent(out) :: Status
          character*3                :: MemOrd
        
          call ZanhuaLowerCase(MemoryOrder,MemOrd)
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
      END SUBROUTINE ZanhuaGetDim

      SUBROUTINE ZanhuaGetSize(MemoryOrder, Size, ep1, ep2, ep3, sp1, sp2, sp3)
          character*(*) ,intent(in)  :: MemoryOrder
          integer, intent(out) :: Size
          integer, intent(in) :: ep1, ep2, ep3, sp1, sp2, sp3
          INTEGER :: ndim, ierr
          ! CHARACTER*80 memord
          ! memord = MemoryOrder
      
          CALL ZanhuaGetDim(TRIM(MemoryOrder), ndim, ierr)
      
          Size = 0
          IF (ndim .EQ. 0) THEN    
            Size = 1
          ELSE IF (ndim .EQ. 1) THEN
            if (ep1 - sp1 > 0) &
            Size = (ep1 - sp1 + 1)
          ELSE IF (ndim .EQ. 2) THEN
            if (ep1 - sp1 > 0 .AND. ep2 - sp2 > 0) &
            Size = (ep1 - sp1 + 1) * (ep2 - sp2 + 1) 
          ELSE IF (ndim .EQ. 3) THEN
            if (ep1 - sp1 > 0 .AND. ep2 - sp2 > 0 .AND. ep3 - sp3 > 0) &
            Size = (ep1 - sp1 + 1) * (ep2 - sp2 + 1) * (ep3 - sp3 + 1)
          ENDIF
          return
      END SUBROUTINE ZanhuaGetSize

      SUBROUTINE ZanhuaGetSizeOfType(Type, Size)
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
      END SUBROUTINE ZanhuaGetSizeOfType

      subroutine BputGetBufferSize(p,switch,grid, zanhua_total_size)
          USE module_domain_type, ONLY : fieldlist
          USE module_domain
          USE module_io
          TYPE(domain), intent(in) :: grid
          TYPE( fieldlist ), intent(inout), POINTER :: p
          INTEGER, INTENT(IN) :: switch
          integer, intent(out):: zanhua_total_size
          integer :: zanhua_grid_size, zanhua_type_size
          INTEGER newswitch, itrace
          
          newswitch = switch
          zanhua_total_size = 0
          DO WHILE ( ASSOCIATED( p ) )
            IF ( p%ProcOrient .NE. 'X' .AND. p%ProcOrient .NE. 'Y' ) THEN
              zanhua_grid_size = 0
              IF (p%Ndim .EQ. 0)  THEN
                ! TODO: check 0 or 1
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newswitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    zanhua_grid_size = 1
                    ! zanhua_num_calls = zanhua_num_calls + 1
                    
                    ! helper; print info
                    zanhua_type_size = 4
                    CALL ZanhuaGetSizeOfType(p%Type, zanhua_type_size)
                    
                    zanhua_total_size = zanhua_total_size + zanhua_grid_size * zanhua_type_size
                    ! CALL ZanhuaGetSize(p%MemoryOrder, zanhua_grid_size, p%ep1, p%ep2, p%ep3, p%sp1, p%sp2, p%sp3)
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 1) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newswitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      ! zanhua_num_calls = zanhua_num_calls + 1
                      CALL ZanhuaGetSize(p%MemoryOrder, zanhua_grid_size, p%ep1, p%ep2, p%ep3, p%sp1, p%sp2, p%sp3)

                      ! helper; print info
                      zanhua_type_size = 4
                      CALL ZanhuaGetSizeOfType(p%Type, zanhua_type_size)
                      ! WRITE(wrf_err_message,*) "**** zanhua_size_of_each_call ", zanhua_grid_size * zanhua_type_size
                      ! CALL wrf_debug( 0 , wrf_err_message )
                      zanhua_total_size = zanhua_total_size + zanhua_grid_size * zanhua_type_size
                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 2) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newswitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) .AND.  &
                    ( .NOT. p%subgrid_x .OR. (p%subgrid_x .AND. grid%sr_x .GT. 0) ) .AND. &
                    ( .NOT. p%subgrid_y .OR. (p%subgrid_y .AND. grid%sr_y .GT. 0) )       &
                  ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      ! zanhua_num_calls = zanhua_num_calls + 1
                      CALL ZanhuaGetSize(p%MemoryOrder, zanhua_grid_size, p%ep1, p%ep2, p%ep3, p%sp1, p%sp2, p%sp3)

                      ! helper; print info
                      zanhua_type_size = 4
                      CALL ZanhuaGetSizeOfType(p%Type, zanhua_type_size)
                      ! WRITE(wrf_err_message,*) "**** zanhua_size_of_each_call ", zanhua_grid_size * zanhua_type_size
                      ! CALL wrf_debug( 0 , wrf_err_message )
                      zanhua_total_size = zanhua_total_size + zanhua_grid_size * zanhua_type_size

                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 3) THEN
                IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams,newswitch)) THEN
                  IF ( in_use_for_config(grid%id,TRIM(p%VarName)) .AND.  &
                    ( .NOT. p%subgrid_x .OR. (p%subgrid_x .AND. grid%sr_x .GT. 0) ) .AND. &
                    ( .NOT. p%subgrid_y .OR. (p%subgrid_y .AND. grid%sr_y .GT. 0) )       &
                  ) THEN
                    IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                      ! zanhua_num_calls = zanhua_num_calls + 1
                      CALL ZanhuaGetSize(p%MemoryOrder, zanhua_grid_size, p%ep1, p%ep2, p%ep3, p%sp1, p%sp2, p%sp3)

                      ! helper; print info
                      zanhua_type_size = 4
                      CALL ZanhuaGetSizeOfType(p%Type, zanhua_type_size)
                      ! WRITE(wrf_err_message,*) "**** zanhua_size_of_each_call ", zanhua_grid_size * zanhua_type_size
                      ! CALL wrf_debug( 0 , wrf_err_message )
                      zanhua_total_size = zanhua_total_size + zanhua_grid_size * zanhua_type_size
                    ENDIF
                  ENDIF
                ENDIF
              ELSE IF (p%Ndim .EQ. 4 .AND. p%scalar_array ) THEN
                IF (switch.EQ.restart_only.OR.p%Ntl/100.EQ.mod(p%Ntl,100)) THEN
                  DO itrace = PARAM_FIRST_SCALAR , p%num_table(grid%id)
                    
                    IF ((p%Restart.AND.switch.EQ.restart_only).OR.on_stream(p%streams_table(grid%id,itrace)%stream,newswitch)) THEN
                      ! zanhua_num_calls = zanhua_num_calls + 1
                      CALL ZanhuaGetSize(p%MemoryOrder, zanhua_grid_size, p%ep1, p%ep2, p%ep3, p%sp1, p%sp2, p%sp3)

                      ! helper; print info
                      zanhua_type_size = 4
                      CALL ZanhuaGetSizeOfType(p%Type, zanhua_type_size)
                      ! WRITE(wrf_err_message,*) "**** zanhua_size_of_each_call ", zanhua_grid_size * zanhua_type_size
                      ! CALL wrf_debug( 0 , wrf_err_message )
                      zanhua_total_size = zanhua_total_size + zanhua_grid_size * zanhua_type_size
                    ENDIF
                  ENDDO
                  
                ENDIF
              ENDIF

            ENDIF
            p => p%next
          ENDDO
        
      end subroutine BputGetBufferSize

  end module