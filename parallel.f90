!>
!! @brief module for parallel processing.
!! @file  parallel.f90
!!

!>
!! @brief module for parallel processing.
module parallel
   use mpi
   implicit none

   !! variables for MPI process.
   integer :: mpi_size !< total number of mpi processes.
   integer :: mpi_rank !< rank of this mpi processes.
   integer :: mpi_ierr !< infomation from MPI subroutines.
   logical :: mpi_parent !< parent process or not.

   !! variables for MPI process grid.
   integer, parameter :: mpi_dim=3
   integer :: mpi_comm_cart            !< communicator among cartesian processes grid.
   integer :: mpi_comm_line(mpi_dim)   !< communicator among processes along each axis.
   integer :: mpi_size_cart(mpi_dim)   !< sizes of cartesian processes grid.
   integer :: mpi_rank_cart(mpi_dim)   !< ranks in cartesian processes grid.
   integer :: mpi_halo_send(2,mpi_dim) !< shape of each halo grid for sending.
   integer :: mpi_halo_recv(2,mpi_dim) !< shape of each halo grid for recving.
   integer :: mpi_body                 !< shape of body grid
   integer :: mpi_body_g               !< shape of local body grid in global grid.

contains
   !>
   !! @brief init MPI and set some variables for MPI process.
   subroutine MPI__Initialize
      call MPI_Init( mpi_ierr )

      call MPI_Comm_size( MPI_COMM_WORLD, mpi_size, mpi_ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, mpi_rank, mpi_ierr )

      mpi_parent = mpi_rank == 0
   end subroutine MPI__Initialize

   !>
   !! @brief divide grid and setup MPI process grid.
   !! @param [out] nx    number of local grid points along x-axis in this process.
   !! @param [out] ny    number of local grid points along y-axis in this process.
   !! @param [out] nz    number of local grid points along z-axis in this process.
   !! @param [out] ex    edges along x-axis of local grid are on the edge of global grid or not.
   !! @param [out] ey    edges along y-axis of local grid are on the edge of global grid or not.
   !! @param [out] ez    edges along z-axis of local grid are on the edge of global grid or not.
   !! @param [out] sx    starting global grid point along x-axis in this process.
   !! @param [out] sy    starting global grid point along y-axis in this process.
   !! @param [out] sz    starting global grid point along z-axis in this process.
   !! @param [in]  nx_g  number of global grid points along x-axis.
   !! @param [in]  ny_g  number of global grid points along y-axis.
   !! @param [in]  nz_g  number of global grid points along z-axis.
   !! @param [in]  nh    number of halo grid points.
   subroutine MPI__DivideGrid( nx, ny, nz, ex, ey, ez, sx, sy, sz, nx_g, ny_g, nz_g, nh )
      integer, intent(out) :: nx, ny, nz
      logical, intent(out) :: ex(2), ey(2), ez(2)
      integer, intent(out) :: sx, sy, sz
      integer, intent(in)  :: nx_g, ny_g, nz_g, nh

      integer :: sxp, syp, szp, nxp, nyp, nzp
      integer :: sxm, sym, szm, nxm, nym, nzm

      !! varibales for finding shape of process grid.
      integer :: size_x, size_y, size_z, size_xyz, size_yz
      integer :: cut, sur, cut_best, sur_best

      !! variables for MPI communication among process grid.
      integer :: grid_size(3), halo_size(3)
      integer :: disp_send(3), disp_recv(3)
      logical, parameter :: periodic(mpi_dim) = [ .false., .false., .false. ] ! x-major 3D process grid.
!!$      logical, parameter :: periodic(mpi_dim) = [ .false., .false., .false. ] ! z-major 3D process grid.

      logical, parameter :: axis_line(mpi_dim,mpi_dim) = &
         ! x-major 3D process grid.
         reshape( (/ &
         .false., .false.,  .true., &
         .false.,  .true., .false., &
         .true. , .false., .false. /), (/mpi_dim,mpi_dim/) )
!!$      ! z-major 3D process grid.
!!$         reshape( (/ &
!!$         .true. , .false., .false., &
!!$         .false.,  .true., .false., &
!!$         .false., .false.,  .true.  /), (/mpi_dim,mpi_dim/) )
      logical, parameter :: ordered = .false.
      integer :: axis, fb

      !! find the best shape of process grid.
      cut_best = 0
      sur_best = 0
      size_xyz = mpi_size
      do size_x=1, size_xyz
         if( mod(size_xyz,size_x) /= 0 ) cycle
         size_yz = size_xyz/size_x

         do size_y=1, size_yz
            if( mod(size_yz,size_y) /= 0 ) cycle
            size_z = size_yz/size_y

            nx = int((nx_g-1)/size_x)+1
            ny = int((ny_g-1)/size_y)+1
            nz = int((nz_g-1)/size_z)+1
            if( nx<=2 .or. ny<=2 .or. nz<=2 ) cycle ! avoid too small grid.

            cut = size_x-1 + size_y-1 + size_z-1
            sur = nx*ny + ny*nz + nz*nx

            if( sur_best==0 .or. sur_best>sur .or. (sur_best==sur .and. cut_best>cut) ) then
               sur_best = sur
               cut_best = cut
               mpi_size_cart(:) = [ size_x, size_y, size_z ]
            end if
         end do ! size_y
      end do ! size_x

      if( sur_best == 0 ) then
         if( mpi_parent ) then
            write(*,*) "Error: failed to divide grid."
            write(*,*) "Error: please change the number of MPI processes from:", mpi_size
         end if
         call MPI_Abort( MPI_COMM_WORLD, 1, mpi_ierr )
      end if

      !! setup process grid.
      call MPI_Cart_create( MPI_COMM_WORLD, mpi_dim, &
           [ mpi_size_cart(3), mpi_size_cart(2), mpi_size_cart(1) ], & ! x-major 3D process grid.
           periodic, ordered, mpi_comm_cart, mpi_ierr )
!!$      call MPI_Cart_create( MPI_COMM_WORLD, mpi_dim, &
!!$           mpi_size_cart, & ! z-major 3D process grid.
!!$           periodic, ordered, mpi_comm_cart, mpi_ierr )

!!$      if( mpi_parent ) then
!!$         write(*,"(a,i6)")      "Info: number of MPI processes:    ", mpi_size
!!$         write(*,"(a,3(x,i4))") "Info: shape  of MPI process grid: ", mpi_size_cart(:)
!!$         write(*,"(a,3(x,i4))") "Info: number of global grid point:", nx_g, ny_g, nz_g
!!$      end if

      !! setup MPI communicators among process grids along x, y, or z-axis.
      do axis=1, mpi_dim
         call MPI_Cart_sub ( mpi_comm_cart, axis_line(:,axis), mpi_comm_line(axis), mpi_ierr )
         call MPI_Comm_rank( mpi_comm_line(axis), mpi_rank_cart(axis), mpi_ierr )
         call MPI_Comm_size( mpi_comm_line(axis), mpi_size_cart(axis), mpi_ierr )
      end do

      !! setup local grid in this process.
      nx = int(nx_g/mpi_size_cart(1)) + merge( 1, 0, mpi_rank_cart(1)<mod(nx_g,mpi_size_cart(1)) )
      ny = int(ny_g/mpi_size_cart(2)) + merge( 1, 0, mpi_rank_cart(2)<mod(ny_g,mpi_size_cart(2)) )
      nz = int(nz_g/mpi_size_cart(3)) + merge( 1, 0, mpi_rank_cart(3)<mod(nz_g,mpi_size_cart(3)) )

      ex = [ mpi_rank_cart(1) == 0, mpi_rank_cart(1) == mpi_size_cart(1)-1 ]
      ey = [ mpi_rank_cart(2) == 0, mpi_rank_cart(2) == mpi_size_cart(2)-1 ]
      ez = [ mpi_rank_cart(3) == 0, mpi_rank_cart(3) == mpi_size_cart(3)-1 ]

      sx = int(nx_g/mpi_size_cart(1))*mpi_rank_cart(1) + min( mpi_rank_cart(1), mod(nx_g,mpi_size_cart(1)) )
      sy = int(ny_g/mpi_size_cart(2))*mpi_rank_cart(2) + min( mpi_rank_cart(2), mod(ny_g,mpi_size_cart(2)) )
      sz = int(nz_g/mpi_size_cart(3))*mpi_rank_cart(3) + min( mpi_rank_cart(3), mod(nz_g,mpi_size_cart(3)) )

      if( mpi_rank_cart(1) == 0 ) then
         nxm = 0
         sxm = 0
      else
         nxm = int(nx_g/mpi_size_cart(1)) + merge( 1, 0, mpi_rank_cart(1)-1<mod(nx_g,mpi_size_cart(1)) )
         sxm = sx - nxm
      end if

      if( mpi_rank_cart(1) == mpi_size_cart(1)-1 ) then
         nxp = 0
         sxp = 0
      else
         nxp = int(nx_g/mpi_size_cart(1)) + merge( 1, 0, mpi_rank_cart(1)+1<mod(nx_g,mpi_size_cart(1)) )
         sxp = sx + nx
      end if

      if( mpi_rank_cart(2) == 0 ) then
         nym = 0
         sym = 0
      else
         nym = int(ny_g/mpi_size_cart(2)) + merge( 1, 0, mpi_rank_cart(2)-1<mod(ny_g,mpi_size_cart(2)) )
         sym = sy - nym
      end if

      if( mpi_rank_cart(2) == mpi_size_cart(2)-1 ) then
         nyp = 0
         syp = 0
      else
         nyp = int(ny_g/mpi_size_cart(2)) + merge( 1, 0, mpi_rank_cart(2)+1<mod(ny_g,mpi_size_cart(2)) )
         syp = sy + ny
      end if

      if( mpi_rank_cart(3) == 0 ) then
         nzm = 0
         szm = 0
      else
         nzm = int(nz_g/mpi_size_cart(3)) + merge( 1, 0, mpi_rank_cart(3)-1<mod(nz_g,mpi_size_cart(3)) )
         szm = sz - nzm
      end if

      if( mpi_rank_cart(3) == mpi_size_cart(3)-1 ) then
         nzp = 0
         szp = 0
      else
         nzp = int(nz_g/mpi_size_cart(3)) + merge( 1, 0, mpi_rank_cart(3)+1<mod(nz_g,mpi_size_cart(3)) )
         szp = sz + nz
      end if

!!$      write(*,"('Info: rank: ',i6,'  cart:',3(x,i2),'  grid: ',3(x,i4),'  shift: ',3(x,i4),'  edge:',6(x,l1))") &
!!$         mpi_rank, mpi_rank_cart(:), nx, ny, nz, sx, sy, sz, ex, ey, ez


      !! setup MPI types for communication of halo grid.
      do axis=1, mpi_dim
         do fb=1, 2
            SELECT CASE( axis )
            CASE(1) ! x
               grid_size = [ nh+nx+nh, nh+ny+nh, nh+nz+nh ]
               halo_size = [ nh, nh+ny+nh, nh+nz+nh ]
               if( fb == 1 ) then ! x-forward
                  disp_send = [ nh,    0, 0 ] ! send from here
                  disp_recv = [  0,    0, 0 ] ! recv to here
               else               ! x-backward
                  disp_send = [ nx,    0, 0 ] ! send from here
                  disp_recv = [ nh+nx, 0, 0 ] ! recv to here
               end if
            CASE(2) ! y
               grid_size = [ nh+nx+nh, nh+ny+nh, nh+nz+nh ]
               halo_size = [ nh+nx+nh, nh, nh+nz+nh ]
               if( fb == 1 ) then ! y-forward
                  disp_send = [ 0, nh,    0 ] ! send from here
                  disp_recv = [ 0,  0,    0 ] ! recv to here
               else               ! y-backward
                  disp_send = [ 0, ny,    0 ] ! send from here
                  disp_recv = [ 0, nh+ny, 0 ] ! recv to here
               end if
            CASE(3) ! z
               grid_size = [ nh+nx+nh, nh+ny+nh, nh+nz+nh ]
               halo_size = [ nh+nx+nh, nh+ny+nh, nh ]
               if( fb == 1 ) then ! z-forward
                  disp_send = [ 0, 0, nh    ] ! send from here
                  disp_recv = [ 0, 0,  0    ] ! recv to here
               else               ! z-backward
                  disp_send = [ 0, 0, nz    ] ! send from here
                  disp_recv = [ 0, 0, nh+nz ] ! recv to here
               end if
            end SELECT

            call MPI_Type_create_subarray &
               ( 3, grid_size, halo_size, disp_send, &
               MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
               mpi_halo_send(fb,axis), mpi_ierr )
            call MPI_Type_commit( mpi_halo_send(fb,axis), mpi_ierr )

            call MPI_Type_create_subarray &
               ( 3, grid_size, halo_size, disp_recv, &
               MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
               mpi_halo_recv(fb,axis), mpi_ierr )
            call MPI_Type_commit( mpi_halo_recv(fb,axis), mpi_ierr )

         end do ! fb
      end do ! axis

      !! setup MPI types for communication of body grid in global grid.
      call MPI_Type_create_subarray( 3, &
         [nx_g,ny_g,nz_g], [nx,ny,nz], [sx,sy,sz], &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_body_g, mpi_ierr )
      call MPI_Type_commit( mpi_body_g, mpi_ierr )

      !! setup MPI types for communication of body grid in local grid.
      call MPI_Type_create_subarray( 3, &
         grid_size, [nx,ny,nz], [nh,nh,nh], &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_body, mpi_ierr )
      call MPI_Type_commit( mpi_body, mpi_ierr )
   end subroutine MPI__DivideGrid

   !>
   !! @brief communicate halo grid among process grid.
   !! @param [in,out] field real scaler field whose halo grid is communicated.
   subroutine MPI__SendrecvHalo( field )
      real(8), intent(inout) :: field(*)

      integer, parameter :: leng(2) = [ 1, 1 ]
      integer(MPI_ADDRESS_KIND), parameter :: disp(2) = [ 0, 0 ] 
      integer :: axis

      do axis=1, mpi_dim
         if( mpi_size_cart(axis) == 1 ) cycle

         call MPI_Neighbor_alltoallw &
            ( field, leng, disp, mpi_halo_send(:,axis), &
            field, leng, disp, mpi_halo_recv(:,axis), &
            mpi_comm_line(axis), mpi_ierr )
      end do
   end subroutine MPI__SendrecvHalo

   !>
   !! @brief read field on body grid from a global file.
   !! @param [in] iunit file unit.
   !! @param [in] offset write position.
   !! @param [out] field real scaler field whose body grid is read from a global file.
   subroutine MPI__File_read_body( iunit, offset, field )
      integer, intent(in) :: iunit
      integer(8), intent(in) :: offset
      real(8), intent(out) :: field(*)

      call MPI_File_set_view( iunit, offset, MPI_DOUBLE_PRECISION, mpi_body_g, "native", MPI_INFO_NULL, mpi_ierr )
      call MPI_File_read_all( iunit, field, 1, mpi_body, MPI_STATUS_IGNORE, mpi_ierr )
   end subroutine MPI__File_read_body

   !>
   !! @brief write field on body grid to a global file.
   !! @param [in] ounit file unit.
   !! @param [in] offset read position.
   !! @param [in] field real scaler field whose body grid is written to a global file.
   subroutine MPI__File_write_body( ounit, offset, field )
      integer, intent(in) :: ounit
      integer(8), intent(in) :: offset
      real(8), intent(in) :: field(*)

      call MPI_File_set_view( ounit, offset, MPI_DOUBLE_PRECISION, mpi_body_g, "native", MPI_INFO_NULL, mpi_ierr )
      call MPI_File_write_all( ounit, field, 1, mpi_body, MPI_STATUS_IGNORE, mpi_ierr )
   end subroutine MPI__File_write_body

   !>
   !! @brief finalize MPI and clear some variables for MPI process grid.
   subroutine MPI__Finalize
      integer :: axis, fb

      call MPI_Type_free( mpi_body, mpi_ierr )
      call MPI_Type_free( mpi_body_g, mpi_ierr )

      do axis=1, mpi_dim
         do fb=1, 2
            call MPI_Type_free( mpi_halo_send(fb,axis), mpi_ierr )
            call MPI_Type_free( mpi_halo_recv(fb,axis), mpi_ierr )
         end do
         call MPI_Comm_free( mpi_comm_line(axis), mpi_ierr )
      end do

      call MPI_Comm_free( mpi_comm_cart, mpi_ierr )

      call MPI_Finalize( mpi_ierr )
   end subroutine MPI__Finalize

end module parallel
