!>
!! @brief Phase-field acid wormholing model (with restart option)
!! @file  wPFAW3D.f90
!! @date  Last updated on 2022.01.26
!!

!> @mainpage
!!
!! - wPFAW3D.f90 : wPFAW3D program.
!! - parallel.f90 : parallel module.
!!

!! @cond no_Doxygen
#include <petsc/finclude/petscsnes.h>
!! @endcond

!>
!! @brief module of wPFAW3D calculation.
module wPFAW3D
   use parallel
   use petscsnes
   use petscksp
   implicit none

   !>
   !! @brief structure of calculation parameters.
   type Param_type
      integer :: nx_g       !< number of global grid points along x-axis.
      integer :: ny_g       !< number of global grid points along y-axis.
      integer :: nz_g       !< number of global grid points along z-axis.
      integer :: ngrid_g    !< total  of global grid points.
      integer :: nx         !< number of local grid points along x-axis in this process.
      integer :: ny         !< number of local grid points along y-axis in this process.
      integer :: nz         !< number of local grid points along z-axis in this process.
      integer :: ngrid      !< total  of local grid points.
      logical :: ex(2)      !< edges along x-axis of local grid are on the edge of global grid or not.
      logical :: ey(2)      !< edges along y-axis of local grid are on the edge of global grid or not.
      logical :: ez(2)      !< edges along z-axis of local grid are on the edge of global grid or not.
      integer :: sx         !< starting global grid point along x-axis in this process.
      integer :: sy         !< starting global grid point along y-axis in this process.
      integer :: sz         !< starting global grid point along z-axis in this process.
      real(8) :: dtime      !< time step [sec].
      real(8) :: tfac       !< factor for dtime.
      real(8) :: time_start !< time at start.
      real(8) :: time_end   !< time at end.
      real(8) :: time_total !< total time from start to end.
      integer :: n_res_out  !< output frequency of result.
      integer :: n_dat_out  !< output frequency of progress.
      real(8) :: dtime_res  !< output interval of result.
      real(8) :: dtime_dat  !< output interval of progress.
      real(8) :: epsilon    !< interface width [m].
      real(8) :: dx         !< grid spacing [m].
      real(8) :: lx         !< x length of grid [m].
      real(8) :: ly         !< y length of grid [m].
      real(8) :: lz         !< z length of grid [m].
      real(8) :: lambda     !< coupling between phi and c [-]
      real(8) :: porosity   !< porosity of rock [-].
      real(8) :: rhos       !< density of rock [kg/m3].
      real(8) :: rhoa       !< density of acid [kg/m3].
      real(8) :: viscosity  !< viscosity of acid [Pa*s].
      real(8) :: cwtp       !< acid concentration [kg HCl/kg solution].
      real(8) :: diffusion  !< diffusion coefficient [m2/s].
      real(8) :: b100       !< gravimetric dissolving power [kg HCl/kg solution].
      real(8) :: ral        !< longitudinal dispersivity [-].
      real(8) :: rat        !< transverse dispersivity [-].
      real(8) :: permeab_hole !< permeability of wormholed rock [m2]
      real(8) :: alpha      !< acid spend rate in dissolution [-].
      real(8) :: pini       !< initial pressure [Pa].
      integer :: ix_ref     !< reference grid x-index.
      integer :: iy_ref     !< reference grid y-index.
      integer :: iz_ref     !< reference grid z-index.
      integer :: nwell      !< total number of well conditions.
      character(len=256) :: filename_well !< file name of input well data.
      character(len=256) :: filename_init !< file name of initial field data to input.
      character(len=256) :: filename_last !< file name of last field data to output.
   end type Param_type

   type(Param_type) :: param !< variable of parameter structure.

   !>
   !! @brief structure of well conditions at a well grid point.
   type Well_type
      integer :: ix_g       !< global grid x-index of this well grid point.
      integer :: iy_g       !< global grid y-index of this well grid point.
      integer :: iz_g       !< global grid z-index of this well grid point.
      real(8) :: pressure   !< fixed pressure at boundary well grid point [Pa].
      real(8) :: injection  !< injection at inlet well grid point [dVol/s].
      real(8) :: phase      !< fixed phase at boundary or inlet well grid point [-].
      real(8) :: concent    !< fixed conentration at boundary or inlet well grid point [-].
   end type Well_type

   type(Well_type), allocatable :: well(:) !< array of all well conditions.

   !! arrays of fields.
   integer, allocatable :: lw(:,:,:)           !< index of well condition at each grid points.
   real(8), allocatable :: permeab_rock(:,:,:) !< permeability of original rock [m2].
   real(8), allocatable :: permeab(:,:,:,:)    !< permeability of rock [m2]. 0:grid point, 1:x-mid, 2:y-mid, 3:z-mid.
   real(8), allocatable :: pressure(:,:,:)     !< pressure [Pa].
   real(8), allocatable :: velocity(:,:,:,:)   !< velocity of acid [m/s]. 1: vx at x-mid, 2: vy at y-mid, 3:vz at z-mid.
   real(8), allocatable :: phase(:,:,:)        !< phase-field variable [-]. -1:solid, +1:liquid.
   real(8), allocatable :: concent(:,:,:)      !< concentration of acid [-].
   real(8), allocatable :: dphasedt(:,:,:)     !< time derivative of phase-field variable [1/s].
   real(8), allocatable :: dconcentdt(:,:,:)   !< time derivative of concentration of acid [1/s].

contains

   !>
   !! @brief load input file.
   !! @param [in] filename file name of input file.
   subroutine loadInput( filename )
      !--- Dummy arguments.
      character(len=*),intent(in) :: filename
      !--- Local variables.
      integer :: iunit, status
      real(8) :: b, kc, wval_, dxw
      integer :: ix, iy, iz
      integer :: iw, wpq_
      real(8) :: init_permeab, init_phase, init_concent

      init_permeab =  0.0d0
      init_phase   = -1.0d0
      init_concent =  0.0d0

      open( newunit=iunit, file=filename, iostat=status, status="old" )
      if( status/=0 ) then
         write(*,*) "Error: can not open input file: ", trim(filename)
         call exit(1)
      end if

      !--- Body of input
      ! Section(1): General data
      read(iunit,*); read(iunit,*); read(iunit,*); read(iunit,*)
      read(iunit,*); read(iunit,*) param%nx_g, param%ny_g, param%nz_g
      read(iunit,*); read(iunit,*) param%dtime, param%tfac
      read(iunit,*); read(iunit,*) param%time_start, param%time_end, param%n_res_out, param%n_dat_out

      ! Section(2): Phase-field parameters
      read(iunit,*)
      read(iunit,*); read(iunit,*) param%epsilon, dxw
      read(iunit,*); read(iunit,*) param%lambda

      ! Section(3): Rock and fluids data
      read(iunit,*)
      read(iunit,*); read(iunit,*) param%porosity
      read(iunit,*); read(iunit,*) param%rhos, param%rhoa
      read(iunit,*); read(iunit,*) param%viscosity
      read(iunit,*); read(iunit,*) param%cwtp
      read(iunit,*); read(iunit,*) param%diffusion

      ! Section(4): Other parameters
      read(iunit,*)
      read(iunit,*); read(iunit,*) param%b100
      read(iunit,*); read(iunit,*) param%ral, param%rat
      read(iunit,*); read(iunit,*) param%permeab_hole
      read(iunit,*); read(iunit,*) param%pini
      read(iunit,*); read(iunit,*) param%ix_ref, param%iy_ref, param%iz_ref

      ! set related parameters.
      param%time_total = param%time_end - param%time_start
      param%dtime_res = param%time_total/param%n_res_out
      param%dtime_dat = param%time_total/param%n_dat_out
      param%dx = dxw*param%epsilon
      param%lx = param%dx*param%nx_g
      param%ly = param%dx*param%ny_g
      param%lz = param%dx*param%nz_g
      b  = param%rhoa/param%rhos
      kc = param%cwtp*param%b100 / (1.0d0-param%porosity)
      param%alpha = 0.5d0/b/kc

      ! Section(5): Well data
      read(iunit,*)
      read(iunit,*); read(iunit,*) param%nwell

      ! set well conditions at well grid points.
      allocate( well(param%nwell) )
      do iw=1, param%nwell
         read(iunit,*) well(iw)%ix_g, well(iw)%iy_g, well(iw)%iz_g, &
            wpq_, wval_, well(iw)%phase, well(iw)%concent

         select case( wpq_ ) !wpqが1か-1で場合分け
         case(+1) ! fixed pressure at boundary well grid point.
            well(iw)%pressure  = wval_
            well(iw)%injection = 0.0d0
         case(-1) ! injection at inlet well grid point.
            well(iw)%pressure  = 0.0d0
            well(iw)%injection = wval_
         case default
            write(*,*) "Error: invalid well data in input file: ", trim(filename)
            call exit(1)
         end select
      end do

      ! Section(6): Initial conditions
      read(iunit,*); read(iunit,*)
      read(iunit,*) init_permeab
      read(iunit,*)
      read(iunit,*) init_phase
      read(iunit,*)
      read(iunit,*) init_concent

      close(iunit)      

      call MPI__DivideGrid( param%nx, param%ny, param%nz, param%ex, param%ey, param%ez, &
         param%sx, param%sy, param%sz, param%nx_g, param%ny_g, param%nz_g, 1 )

      param%ngrid    = param%nx  *param%ny  *param%nz
      param%ngrid_g  = param%nx_g*param%ny_g*param%nz_g

      ! allocate fields.
      allocate( lw(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( permeab_rock(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( permeab(0:3,0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( pressure(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( velocity(3,0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( phase(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( concent(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( dphasedt(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( dconcentdt(0:param%nx+1,0:param%ny+1,0:param%nz+1) )

      ! set well index table.
      lw(:,:,:) = 0
      do iw=1, param%nwell    !グローバル座標からローカル座標ix,iy,izに変換
         ix = well(iw)%ix_g - param%sx
         iy = well(iw)%iy_g - param%sy
         iz = well(iw)%iz_g - param%sz

         if( 0<=ix .and. ix<=param%nx+1 .and. 0<=iy .and. iy<=param%ny+1 .and. 0<=iz .and. iz<=param%nz+1 ) then
            lw(ix,iy,iz) = iw ! index of well grid point in array well.         
         end if
      end do

      ! initialize fields.
      permeab_rock(:,:,:) = init_permeab
      permeab(:,:,:,:)    = 0.0d0
      pressure(:,:,:)     = param%pini
      velocity(:,:,:,:)   = 0.0d0
      phase(:,:,:)        = init_phase
      concent(:,:,:)      = init_concent
      dphasedt(:,:,:)     = 0.0d0
      dconcentdt(:,:,:)   = 0.0d0

      ! set fixed values for phase and concent at well grid points.
      !phase concentに初期条件代入
      do iw=1, param%nwell  
         ix = well(iw)%ix_g - param%sx
         iy = well(iw)%iy_g - param%sy
         iz = well(iw)%iz_g - param%sz

         if( 1<=ix .and. ix<=param%nx .and. 1<=iy .and. iy<=param%ny .and. 1<=iz .and. iz<=param%nz ) then
            phase(ix,iy,iz)   = well(iw)%phase
            concent(ix,iy,iz) = well(iw)%concent
         end if
      end do

      !各行列のひとつ前の値を袖(+1の部分1)で拡大
      call updateHalo( permeab_rock )
      call updateHalo( phase )
      call updateHalo( concent )
   end subroutine loadInput

   !>
   !! @brief load configuration file.
   !! @param [in] filename file name of configuration file.
   subroutine loadConfig( filename )
      !--- Dummy arguments.
      character(len=*), intent(in) :: filename
      !--- Local variables.
      integer :: iunit, status
      character(len=1024) :: buf, key
      real(8) :: b, kc, dxw, wval_
      integer :: i, iw, wpq_
      integer :: ix, iy, iz
      real(8) :: init_permeab, init_phase, init_concent

      init_permeab =  0.0d0
      init_phase   = -1.0d0
      init_concent =  0.0d0
      param%filename_init = ""

      open( newunit=iunit, file=filename, iostat=status, status="old" )
      if( status/=0 ) then
         write(*,*) "Error: can not open config file: ", trim(filename)
         call exit(1)
      end if

      param%dtime = 0.1000d-06

      do
         read(iunit,'(a)',iostat=status) buf
         if( status/=0 ) exit ! end of file.

         if( buf(1:1)=="#" ) cycle ! skip comment line

         i = index(buf,"#")
         if( i>0 ) buf = buf(1:i-1) ! cut comment words.

         if( len_trim(buf)==0 ) cycle ! skip empty line.

         i = index(buf,"=")
         if( i>0 ) buf(i:i) = " " ! replace "=" by a white space.

         read(buf,*,iostat=status) key ! read the first word as keyword.
         if( status/=0 ) goto 1000

         select case( key )
         case("grid")
            read(buf,*,iostat=status) key, param%nx_g, param%ny_g, param%nz_g
            if( status/=0 ) goto 1000
         case("dt")
            read(buf,*,iostat=status) key, param%dtime
            if( status/=0 ) goto 1000
         case("tfac")
            read(buf,*,iostat=status) key, param%tfac
            if( status/=0 ) goto 1000
         case("time_total")
            read(buf,*,iostat=status) key, param%time_total
            if( status/=0 ) goto 1000
         case("gid_output")
            read(buf,*,iostat=status) key, param%n_res_out
            if( status/=0 ) goto 1000
         case("scr_output")
            read(buf,*,iostat=status) key, param%n_dat_out
            if( status/=0 ) goto 1000
         case("w")
            read(buf,*,iostat=status) key, param%epsilon
            if( status/=0 ) goto 1000
         case("dxw")
            read(buf,*,iostat=status) key, dxw
            if( status/=0 ) goto 1000
         case("lambda")
            read(buf,*,iostat=status) key, param%lambda
            if( status/=0 ) goto 1000
         case("poro")
            read(buf,*,iostat=status) key, param%porosity
            if( status/=0 ) goto 1000
         case("rhos")
            read(buf,*,iostat=status) key, param%rhos
            if( status/=0 ) goto 1000
         case("rhoa")
            read(buf,*,iostat=status) key, param%rhoa
            if( status/=0 ) goto 1000
         case("visco")
            read(buf,*,iostat=status) key, param%viscosity
            if( status/=0 ) goto 1000
         case("cwtp")
            read(buf,*,iostat=status) key, param%cwtp
            if( status/=0 ) goto 1000
         case("D")
            read(buf,*,iostat=status) key, param%diffusion
            if( status/=0 ) goto 1000
         case("b100")
            read(buf,*,iostat=status) key, param%b100
            if( status/=0 ) goto 1000
         case("ral")
            read(buf,*,iostat=status) key, param%ral
            if( status/=0 ) goto 1000
         case("rat")
            read(buf,*,iostat=status) key, param%rat
            if( status/=0 ) goto 1000
         case("kwh")
            read(buf,*,iostat=status) key, param%permeab_hole
            if( status/=0 ) goto 1000
         case("pini")
            read(buf,*,iostat=status) key, param%pini
            if( status/=0 ) goto 1000
         case("ref")
            read(buf,*,iostat=status) key, param%ix_ref, param%iy_ref, param%iz_ref
            if( status/=0 ) goto 1000
         case("well")
            read(buf,*,iostat=status) key, param%filename_well
            if( status/=0 ) goto 1000
         case("init")
            read(buf,*,iostat=status) key, init_permeab, init_phase, init_concent
            if( status/=0 ) then
               read(buf,*,iostat=status) key, param%filename_init
               if( status/=0 ) goto 1000
            end if
         case("last")
            read(buf,*,iostat=status) key, param%filename_last
            if( status/=0 ) goto 1000
         case default
            write(*,*) "Error: unknown key in config file: ", trim(filename)
            write(*,*) "Error: unknown key: ", trim(key)
            goto 1000
         end select
      end do

      close(iunit)

      ! set related parameters.
      param%dtime_res = param%time_total/param%n_res_out
      param%dtime_dat = param%time_total/param%n_dat_out
      param%dx = dxw*param%epsilon
      param%lx = param%dx*param%nx_g
      param%ly = param%dx*param%ny_g
      param%lz = param%dx*param%nz_g
      b  = param%rhoa/param%rhos
      kc = param%cwtp*param%b100/(1.0d0 - param%porosity)
      param%alpha = 0.5d0/b/kc

      ! read well file.
      open( newunit=iunit, file=param%filename_well, iostat=status, status="old" )
      if( status/=0 ) then
         write(*,*) "Error: can not open well file: ", trim(param%filename_well)
         call exit(1)
      end if

      read(iunit,'(a)',iostat=status) buf
      read(iunit,'(a)',iostat=status) buf
      read(buf,*,iostat=status) param%nwell

      ! set well conditions at well grid points.
      allocate( well(param%nwell) )
      do iw=1, param%nwell
         read(iunit,'(a)',iostat=status) buf
         read(buf,*) well(iw)%ix_g, well(iw)%iy_g, well(iw)%iz_g, &
            wpq_, wval_, well(iw)%phase, well(iw)%concent

         select case( wpq_ )
         case(+1) ! fixed pressure at boundary well grid point.
            well(iw)%pressure  = wval_
            well(iw)%injection = 0.0d0
         case(-1) ! injection at inlet well grid point.
            well(iw)%pressure  = 0.0d0
            well(iw)%injection = wval_
         case default
            write(*,*) "Error: invalid data in well file: ", trim(param%filename_well)
            call exit(1)
         end select
      end do

      close(iunit)

      call MPI__DivideGrid( param%nx, param%ny, param%nz, param%ex, param%ey, param%ez, &
         param%sx, param%sy, param%sz, param%nx_g, param%ny_g, param%nz_g, 1 )

      param%ngrid    = param%nx  *param%ny  *param%nz
      param%ngrid_g  = param%nx_g*param%ny_g*param%nz_g

      ! allocate fields.
      allocate( lw(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( permeab_rock(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( permeab(0:3,0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( pressure(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( velocity(3,0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( phase(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( concent(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( dphasedt(0:param%nx+1,0:param%ny+1,0:param%nz+1) )
      allocate( dconcentdt(0:param%nx+1,0:param%ny+1,0:param%nz+1) )

      ! set well index table.
      lw(:,:,:) = 0
      do iw=1, param%nwell
         ix = well(iw)%ix_g - param%sx
         iy = well(iw)%iy_g - param%sy
         iz = well(iw)%iz_g - param%sz

         if( 0<=ix .and. ix<=param%nx+1 .and. 0<=iy .and. iy<=param%ny+1 .and. 0<=iz .and. iz<=param%nz+1 ) then
            lw(ix,iy,iz) = iw ! index of well grid point in array well.
         end if
      end do

      ! initialize fields.
      permeab_rock(:,:,:) = init_permeab
      permeab(:,:,:,:)    = 0.0d0
      pressure(:,:,:)     = param%pini
      velocity(:,:,:,:)   = 0.0d0
      phase(:,:,:)        = init_phase
      concent(:,:,:)      = init_concent
      dphasedt(:,:,:)     = 0.0d0
      dconcentdt(:,:,:)   = 0.0d0

      ! set fixed values for phase and concent at well grid points.
      do iw=1, param%nwell
         ix = well(iw)%ix_g - param%sx
         iy = well(iw)%iy_g - param%sy
         iz = well(iw)%iz_g - param%sz

         if( 1<=ix .and. ix<=param%nx .and. 1<=iy .and. iy<=param%ny .and. 1<=iz .and. iz<=param%nz ) then
            phase(ix,iy,iz)   = well(iw)%phase
            concent(ix,iy,iz) = well(iw)%concent
         end if
      end do

      if( param%filename_init /= "" ) then
         call loadField( param%time_start )
      else
         param%time_start = 0.0d0
      end if

      ! reset time end.
      param%time_end = param%time_start + param%time_total

      call updateHalo( permeab_rock )
      call updateHalo( phase )
      call updateHalo( concent )

      return

1000  continue
      close(iunit)
      write(*,*) "Error: broken line in config file: ", trim(filename)
      write(*,*) "Error: broken line: ", trim(buf)
      call exit(1)
   end subroutine loadConfig

   !>
   !! @brief load field file.
   subroutine loadField( time )
      !--- Dummy arguments.
      real(8), intent(out) :: time
      !--- Local variables.
      integer :: iunit, status
      integer(8) :: offset

      type Header_type
         character(len=8) :: magic
         character(len=8) :: version
         real(8)    :: time
         integer(2) :: nx_g, ny_g, nz_g, dumm
      end type Header_type
      type(Header_type) :: head

      call MPI_File_open( MPI_COMM_WORLD, param%filename_init, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, iunit, status )
      if( status/=0 ) then
         write(*,*) "Error: can not open field file: ", trim(param%filename_init)
         call exit(1)
      end if

      offset = 0
      call MPI_File_read_all( iunit, head, int(sizeof(head)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
      offset = offset + sizeof(head)

      if( head%magic /= "wPFAWHFG" ) then
         write(*,*) "Error: invalid magic in field file: ", trim(param%filename_init)
         call exit(1)
      end if

      time = head%time

      if( head%nx_g /= param%nx_g .or. head%ny_g /= param%ny_g .or. head%nz_g /= param%nz_g ) then
         write(*,*) "Error: invalid grid in field file: ", trim(param%filename_init)
         call exit(1)
      end if

      call MPI__File_read_body( iunit, offset, permeab_rock )
      offset = offset + param%ngrid_g*sizeof(permeab_rock(1,1,1))

      call MPI__File_read_body( iunit, offset, phase )
      offset = offset + param%ngrid_g*sizeof(phase(1,1,1))

      call MPI__File_read_body( iunit, offset, concent )
      offset = offset + param%ngrid_g*sizeof(concent(1,1,1))

      call MPI_File_close( iunit, mpi_ierr )
   end subroutine loadField

   !>
   !! @brief calculate permeability.
   subroutine calcPermeability
      !--- Local variables.
      integer :: ix, iy, iz
      real(8) :: f

      !-- Update permeabilities at grid points.
      do iz=0, param%nz+1
         do iy=0, param%ny+1
            do ix=0, param%nx+1
               f = 0.5d0*(1.0d0+phase(ix,iy,iz))
               permeab(0,ix,iy,iz) = param%permeab_hole*f + permeab_rock(ix,iy,iz)*(1.0d0-f) ! grid point.
            end do
         end do
      end do

      !-- Update permeabilities at middle of grid points.
      do iz=0, param%nz
         do iy=0, param%ny
            do ix=0, param%nx
               permeab(1,ix,iy,iz) = harmonic_mean( permeab(0,ix,iy,iz), permeab(0,ix+1,iy,iz) ) ! x-mid.
               permeab(2,ix,iy,iz) = harmonic_mean( permeab(0,ix,iy,iz), permeab(0,ix,iy+1,iz) ) ! y-mid.
               permeab(3,ix,iy,iz) = harmonic_mean( permeab(0,ix,iy,iz), permeab(0,ix,iy,iz+1) ) ! z-mid.
            end do
         end do
      end do

   contains
      !>
      !! @brief calculate harmonic mean of given two values.
      !! @param [in] permeab1 a value.
      !! @param [in] permeab2 another value.
      !! @retval harmonic mean.
      function harmonic_mean( permeab1, permeab2 ) result(mean)
         real(8), intent(in) :: permeab1, permeab2
         real(8) :: mean

         mean = 2.0d0/(1.0d0/permeab1 + 1.0d0/permeab2)
      end function harmonic_mean
   end subroutine calcPermeability

   !>
   !! @brief solve pressure.
   subroutine solvePressure
      !--- Local variables.
      ! grid indexes.
      integer :: ix, iy, iz
      integer :: ix_g, iy_g, iz_g
      integer :: nnz, n
      ! tables.
      integer, allocatable :: irn(:), jcn(:)
      ! matrix and vectors.
      real(8), allocatable :: mat_A(:), vec_b(:,:,:), vec_x(:,:,:)
      ! matrix elements.
      real(8) :: txp, txm, typ, tym, tzp, tzm, tc
      ! PETSc parameters.
      PetscInt,  parameter :: max_iterations = 1000
      PetscReal, parameter :: tolerance = 1.0d-10
      PetscInt  :: iterations
      PetscReal :: residual
      PetscInt  :: start_scat, end_scat
      PetscInt, allocatable :: indices_send(:,:,:), indices_recv(:,:,:)
      IS  :: is_send, is_recv
      KSP :: ksp  ! linear solver context
      PC  :: pc   ! preconditioner context
      Mat :: Apetsc
      Vec :: bpetsc, xpetsc, xscat
      VecScatter :: scat  ! scatter context

      !--- Count number of non-zero elements in A matrix.
      nnz=0
      do iz=1,param%nz
         iz_g = iz + param%sz
         do iy=1,param%ny
            iy_g = iy + param%sy
            do ix=1,param%nx
               ix_g = ix + param%sx

               if( .false. ) then
               else if( onBoundary(ix,iy,iz) ) then
                  ! count only diagonal element.
                  nnz=nnz+1
               else
                  ! if the Z prev point is inside the global grid.
                  if( iz_g-1>=1 ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix,iy,iz-1) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix,iy,iz-1).
                        nnz=nnz+1
                     end if
                  end if

                  ! if the Y prev point is inside the global grid.
                  if( iy_g-1>=1 ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix,iy-1,iz) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix,iy-1,iz).
                        nnz=nnz+1
                     end if
                  end if

                  ! if the X prev point is inside the global grid.
                  if( ix_g-1>=1 ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix-1,iy,iz) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix-1,iy,iz).
                        nnz=nnz+1
                     end if
                  end if

                  ! count diagonal element between (ix,iy,iz)-(ix,iy,iz).
                  nnz=nnz+1

                  ! if the X next point is inside the global grid.
                  if( ix_g+1<=param%nx_g  ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix+1,iy,iz) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix+1,iy,iz).
                        nnz=nnz+1
                     end if
                  end if

                  ! if the Y next point is inside the global grid.
                  if( iy_g+1<=param%ny_g  ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix,iy+1,iz) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix,iy+1,iz).
                        nnz=nnz+1
                     end if
                  end if

                  ! if the Z next point is inside the global grid.
                  if( iz_g+1<=param%nz_g  ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix,iy,iz+1) ) then
                        ! count offdiagonal element between (ix,iy,iz)-(ix,iy,iz+1).
                        nnz=nnz+1
                     end if
                  end if
               end if

            end do
         end do
      end do

      !--- Allocate index tables.
      allocate( irn(nnz), jcn(nnz) )

      !--- Allocate matrix and vectors.
      allocate( mat_A(nnz), vec_b(param%nx,param%ny,param%nz), vec_x(param%nx,param%ny,param%nz) )

      !--- Construct index tables.
      nnz=0
      do iz=1,param%nz
         iz_g = iz + param%sz
         do iy=1,param%ny
            iy_g = iy + param%sy
            do ix=1,param%nx
               ix_g = ix + param%sx

               if( .false. ) then
               else if( onBoundary(ix,iy,iz) ) then
                  ! set tables for only diagonal element.
                  nnz=nnz+1
                  irn(nnz) = index3D_g(ix, iy, iz )-1
                  jcn(nnz) = index3D_g(ix, iy, iz )-1
               else
                  ! if the Z prev point is inside the global grid.
                  if( iz_g-1>=1    ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix,iy,iz-1) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix,iy,iz-1).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,iy,iz  )-1
                        jcn(nnz) = index3D_g(ix,iy,iz-1)-1
                     end if
                  end if

                  ! if the Y prev point is inside the global grid.
                  if( iy_g-1>=1    ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix,iy-1,iz) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix,iy-1,iz).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,iy  ,iz)-1
                        jcn(nnz) = index3D_g(ix,iy-1,iz)-1
                     end if
                  end if

                  ! if the X prev point is inside the global grid.
                  if( ix_g-1>=1    ) then
                     ! if the prev point is not fixed.
                     if( .not. onBoundary(ix-1,iy,iz) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix-1,iy,iz).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,  iy,iz)-1
                        jcn(nnz) = index3D_g(ix-1,iy,iz)-1
                     end if
                  end if

                  ! set tables for diagonal element between (ix,iy,iz)-(ix,iy,iz).
                  nnz=nnz+1
                  irn(nnz) = index3D_g(ix, iy, iz )-1
                  jcn(nnz) = index3D_g(ix, iy, iz )-1

                  ! if the X next point is inside the global grid.
                  if( ix_g+1<=param%nx_g ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix+1,iy,iz) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix+1,iy,iz).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,  iy,iz)-1
                        jcn(nnz) = index3D_g(ix+1,iy,iz)-1
                     end if
                  end if

                  ! if the Y next point is inside the global grid.
                  if( iy_g+1<=param%ny_g ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix,iy+1,iz) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix,iy+1,iz).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,iy,  iz)-1
                        jcn(nnz) = index3D_g(ix,iy+1,iz)-1
                     end if
                  end if

                  ! if the Z next point is inside the global grid.
                  if( iz_g+1<=param%nz_g ) then
                     ! if the next point is not fixed.
                     if( .not. onBoundary(ix,iy,iz+1) ) then
                        ! set tables for offdiagonal element between (ix,iy,iz)-(ix,iy,iz+1).
                        nnz=nnz+1
                        irn(nnz) = index3D_g(ix,iy,iz  )-1
                        jcn(nnz) = index3D_g(ix,iy,iz+1)-1
                     end if
                  end if
               end if

            end do
         end do
      end do

      !--- Construct right-hand vector.
      do iz=1,param%nz
         iz_g = iz + param%sz
         do iy=1,param%ny
            iy_g = iy + param%sy
            do ix=1,param%nx
               ix_g = ix + param%sx

               if( .false. ) then
               else if( onBoundary(ix,iy,iz) ) then
                  ! give fixed P through b.
                  vec_b(ix,iy,iz) =  valuePressure(ix,iy,iz)
               else if( onInlet(ix,iy,iz) ) then
                  ! injection.
                  vec_b(ix,iy,iz) = -valueInjection(ix,iy,iz) * (param%dx*param%dx/param%dx**3)
               else
                  ! no injection.
                  vec_b(ix,iy,iz) = 0.0d0
               end if

            end do
         end do
      end do

      !--- Construct and symmetrize A matrix and right-hand b vector.
      nnz=0
      do iz=1,param%nz
         iz_g = iz + param%sz
         do iy=1,param%ny
            iy_g = iy + param%sy
            do ix=1,param%nx
               ix_g = ix + param%sx

               if( .false. ) then
               else if( onBoundary(ix,iy,iz) ) then
                  ! set matrix for only diagonal element.
                  nnz=nnz+1
                  mat_A(nnz) = 1.0d0
               else
                  ! matrix elements.
                  txp = merge( permeab(1,ix,  iy,   iz)/param%viscosity, 0.0d0, ix_g<param%nx_g ) ! (ix+1,iy,iz)
                  txm = merge( permeab(1,ix-1,iy,   iz)/param%viscosity, 0.0d0, ix_g>1  )         ! (ix-1,iy,iz)
                  typ = merge( permeab(2,ix,  iy,   iz)/param%viscosity, 0.0d0, iy_g<param%ny_g ) ! (ix,iy+1,iz)
                  tym = merge( permeab(2,ix,  iy-1, iz)/param%viscosity, 0.0d0, iy_g>1  )         ! (ix,iy-1,iz)
                  tzp = merge( permeab(3,ix,  iy,   iz)/param%viscosity, 0.0d0, iz_g<param%nz_g ) ! (ix,iy,iz+1)
                  tzm = merge( permeab(3,ix,  iy, iz-1)/param%viscosity, 0.0d0, iz_g>1  )         ! (ix,iy,iz-1)
                  tc  = -(txp+txm+typ+tym+tzp+tzm)

                  ! if the Z prev point is inside the global grid.
                  if( iz_g-1>=1 ) then
                     ! if the prev point is fixed.
                     if( onBoundary(ix,iy,iz-1) ) then
                        ! substract row of (ix,iy,iz-1) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix,iy,iz-1)*tzm
                     else
                        ! set matrix for offdiagonal element between (ix,iy,iz)-(ix,iy,iz-1).
                        nnz=nnz+1
                        mat_A(nnz) = tzm
                     end if
                  end if

                  ! if the Y prev point is inside the global grid.
                  if( iy_g-1>=1 ) then
                     ! if the prev point is fixed.
                     if( onBoundary(ix,iy-1,iz) ) then
                        ! substract row of (ix,iy-1,iz) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix,iy-1,iz)*tym
                     else
                        ! set matrix for offdiagonal element between (ix,iy,iz)-(ix,iy-1,iz).
                        nnz=nnz+1
                        mat_A(nnz) = tym
                     end if
                  end if

                  ! if the X prev point is inside the global grid.
                  if( ix_g-1>=1 ) then
                     ! if the prev point is fixed.
                     if( onBoundary(ix-1,iy,iz) ) then
                        ! substract row of (ix-1,iy,iz) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix-1,iy,iz)*txm
                     else
                        ! set matrix for offdiagonal element between (ix,iy)-(ix-1,iy,iz).
                        nnz=nnz+1
                        mat_A(nnz) = txm
                     end if
                  end if

                  ! set matrix for diagonal element between (ix,iy,iz)-(ix,iy,iz).
                  nnz=nnz+1
                  mat_A(nnz) = tc

                  ! if the X next point is inside the global grid.
                  if( ix_g+1<=param%nx_g ) then
                     ! if the next point is fixed.
                     if( onBoundary(ix+1,iy,iz) ) then
                        ! substract row of (ix+1,iy,iz) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix+1,iy,iz)*txp
                     else
                        ! set matrix for offdiagonal element between (ix,iy,iz)-(ix+1,iy,iz).
                        nnz=nnz+1
                        mat_A(nnz) = txp
                     end if
                  end if

                  ! if the Y next point is inside the global grid.
                  if( iy_g+1<=param%ny_g ) then
                     ! if the next point is fixed.
                     if( onBoundary(ix,iy+1,iz) ) then
                        ! substract row of (ix,iy+1,iz) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix,iy+1,iz)*typ
                     else
                        ! set matrix for offdiagonal element between (ix,iy,iz)-(ix,iy+1,iz).
                        nnz=nnz+1
                        mat_A(nnz) = typ
                     end if
                  end if

                  ! if the Z next point is inside the global grid.
                  if( iz_g+1<=param%nz_g ) then
                     ! if the next point is fixed.
                     if( onBoundary(ix,iy,iz+1) ) then
                        ! substract row of (ix,iy,iz+1) from row of (ix,iy,iz).
                        vec_b(ix,iy,iz) = vec_b(ix,iy,iz) - valuePressure(ix,iy,iz+1)*tzp
                     else
                        ! set matrix for offdiagonal element between (ix,iy,iz)-(ix,iy,iz+1).
                        nnz=nnz+1
                        mat_A(nnz) = tzp
                     end if
                  end if
               end if

            end do
         end do
      end do

      !--- Construct PETSc global arrays.
      ! PETSc global A matrix.
      call MatCreate( PETSC_COMM_WORLD, Apetsc, mpi_ierr )
      call MatSetSizes( Apetsc, PETSC_DECIDE, PETSC_DECIDE, param%ngrid_g, param%ngrid_g, mpi_ierr )
      call MatSetFromOptions( Apetsc, mpi_ierr )
      call MatSetUP( Apetsc, mpi_ierr )
      call MatSetOption( Apetsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, mpi_ierr )
      call MatSetOption( Apetsc, MAT_SYMMETRIC, PETSC_TRUE, mpi_ierr )

      do n=1,nnz
         call MatSetValues( Apetsc, 1, irn(n), 1, jcn(n), mat_A(n), INSERT_VALUES, mpi_ierr )
      end do

      call MatAssemblyBegin( Apetsc, MAT_FINAL_ASSEMBLY, mpi_ierr )
      call MatAssemblyEnd( Apetsc, MAT_FINAL_ASSEMBLY, mpi_ierr )

      ! PETSc global b vector.
      call VecCreate( PETSC_COMM_WORLD, bpetsc, mpi_ierr )
      call VecSetSizes( bpetsc, PETSC_DECIDE, param%ngrid_g, mpi_ierr )
      call VecSetFromOptions( bpetsc, mpi_ierr )
      call VecSetUp( bpetsc, mpi_ierr )

      do iz=1,param%nz
         do iy=1,param%ny
            do ix=1,param%nx
               call VecSetValues( bpetsc, 1, index3D_g(ix,iy,iz)-1, vec_b(ix,iy,iz), INSERT_VALUES, mpi_ierr )
            end do
         end do
      end do

      call VecAssemblyBegin( bpetsc, mpi_ierr )
      call VecAssemblyEnd( bpetsc, mpi_ierr )

      ! PETSc global x vector.
      call VecCreate( PETSC_COMM_WORLD, xpetsc, mpi_ierr )
      call VecSetSizes( xpetsc, PETSC_DECIDE, param%ngrid_g, mpi_ierr )
      call VecSetFromOptions( xpetsc, mpi_ierr )
      call VecSetUp( xpetsc, mpi_ierr )

      do iz=1,param%nz
         do iy=1,param%ny
            do ix=1,param%nx
               call VecSetValues( xpetsc, 1, index3D_g(ix,iy,iz)-1, pressure(ix,iy,iz), INSERT_VALUES, mpi_ierr )
            end do
         end do
      end do

      call VecAssemblyBegin( xpetsc, mpi_ierr )
      call VecAssemblyEnd( xpetsc, mpi_ierr )
      
      ! KSP.
      call KSPCreate( PETSC_COMM_WORLD, ksp, mpi_ierr )
      call KSPSetOperators( ksp, Apetsc, Apetsc, mpi_ierr )
      ! precondition.
      call KSPGetPC( ksp, pc, mpi_ierr )
      ! SOR.
      call PCSetType( pc, PCSOR, mpi_ierr )
      call PCSORSetIterations( pc, 10, 1, mpi_ierr )
      call PCSORSetOmega( pc, 1.7d0, mpi_ierr )
      ! KSP
      call KSPSetType( ksp, KSPBCGS, mpi_ierr )
      call KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, max_iterations, mpi_ierr )
      call KSPSetInitialGuessNonzero( ksp, PETSC_TRUE, mpi_ierr )
      call KSPSetFromOptions( ksp, mpi_ierr )

      !--- Solve linear equation A*x=b.
      call KSPSolve( ksp, bpetsc, xpetsc, mpi_ierr )
      if( mpi_ierr/=0 ) stop "ERROR KSPSolve"

!!$      !   output iterations and residual.
!!$      call KSPGetIterationNumber( ksp, iterations, mpi_ierr )
!!$      call KSPGetResidualNorm( ksp, residual, mpi_ierr )
!!$
!!$      if(mpi_rank == 0) then
!!$         write(*,'(a,i5,a,1p,e12.5)') 'BiCGStab  Iterations : ', iterations, ', Residual Norm : ', residual
!!$      end if

      !--- Scatter global x vector to local x vector.
      call VecCreateMPI( PETSC_COMM_WORLD,param%ngrid, PETSC_DECIDE, xscat, mpi_ierr )
      call VecSetFromOptions( xscat, mpi_ierr )
      call VecSet( xscat, 0.0d0, mpi_ierr )
      call VecGetOwnershipRange( xscat, start_scat, end_scat, mpi_ierr )

      allocate( indices_send(param%nx,param%ny,param%nz), indices_recv(param%nx,param%ny,param%nz) )

      do iz=1,param%nz
         do iy=1,param%ny
            do ix=1,param%nx
               indices_send(ix,iy,iz) = index3D_g(ix,iy,iz)-1
               indices_recv(ix,iy,iz) = index3D_l(ix,iy,iz)-1+start_scat
            end do
         end do
      end do
      
      call ISCreateGeneral( PETSC_COMM_WORLD, param%ngrid, indices_send, PETSC_COPY_VALUES, is_send, mpi_ierr )
      call ISCreateGeneral( PETSC_COMM_WORLD, param%ngrid, indices_recv, PETSC_COPY_VALUES, is_recv, mpi_ierr )
      call VecScatterCreate( xpetsc, is_send, xscat, is_recv, scat, mpi_ierr )
      call VecScatterBegin( scat, xpetsc, xscat, INSERT_VALUES, SCATTER_FORWARD, mpi_ierr )
      call VecScatterEnd( scat, xpetsc, xscat, INSERT_VALUES, SCATTER_FORWARD, mpi_ierr )

      do iz=1,param%nz
         do iy=1,param%ny
            do ix=1,param%nx
               call VecGetValues( xscat, 1, indices_recv(ix,iy,iz), vec_x(ix,iy,iz), mpi_ierr )
            end do
         end do
      end do

      deallocate( indices_send, indices_recv )

      pressure(1:param%nx,1:param%ny,1:param%nz) = vec_x(1:param%nx,1:param%ny,1:param%nz)

      call updateHalo( pressure )

      !--- Free PETSc global arrays.
      call MatDestroy( Apetsc, mpi_ierr )
      call VecDestroy( bpetsc, mpi_ierr )
      call VecDestroy( xpetsc, mpi_ierr )
      call VecDestroy( xscat, mpi_ierr )
      call KSPDestroy( ksp, mpi_ierr )
      call ISDestroy( is_send, mpi_ierr )
      call ISDestroy( is_recv, mpi_ierr )
      call VecScatterDestroy( scat, mpi_ierr )

      !--- Free local arrays.
      deallocate( irn, jcn )
      deallocate( mat_A, vec_b, vec_x )
      
   contains
      !>
      !! @brief return serial index of a given grid point in local grid.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along y-axis
      !! @retval serial index in local grid.
      function index3D_l( ix, iy, iz ) result(index)
         integer, intent(in) :: ix, iy, iz
         integer :: index

         index = (ix-1) + (iy-1)*param%nx + (iz-1)*param%nx*param%ny + 1
      end function index3D_l

      !>
      !! @brief return serial index of a given grid point in global grid.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along y-axis
      !! @retval serial index in global grid.
      function index3D_g( ix, iy, iz ) result(index)
         integer, intent(in) :: ix, iy, iz
         integer :: index

         index = (ix+param%sx-1) + (iy+param%sy-1)*param%nx_g + (iz+param%sz-1)*param%nx_g*param%ny_g + 1
      end function index3D_g

      !>
      !! @brief check a given grid point is on boundary or not.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along z-axis
      !! @retval boundary or not.
      function onBoundary( ix, iy, iz ) result(boundary)
         integer, intent(in) :: ix, iy, iz
         integer :: iw
         logical :: boundary

         iw=lw(ix,iy,iz)
         if( iw>0 ) then
            boundary = well(iw)%pressure/=0.0d0
         else
            boundary = .false.
         end if
      end function onBoundary

      !>
      !! @brief return pressure at a given boundary well grid point.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along z-axis
      !! @retval pressure on boundary.
      function valuePressure( ix, iy, iz ) result(pressure)
         integer, intent(in) :: ix, iy, iz
         integer :: iw
         real(8) :: pressure

         iw=lw(ix,iy,iz)
         pressure = well(iw)%pressure
      end function valuePressure

      !>
      !! @brief check a given grid point is on inlet or not.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along z-axis
      !! @retval inlet or not.
      function onInlet( ix, iy, iz ) result(inlet)
         integer, intent(in) :: ix, iy, iz
         integer :: iw
         logical :: inlet

         iw=lw(ix,iy,iz)
         if( iw>0 ) then
            inlet = well(iw)%injection/=0.0d0
         else
            inlet = .false.
         end if
      end function onInlet

      !>
      !! @brief return injection at a given inlet well grid point.
      !! @param [in] ix local index along x-axis
      !! @param [in] iy local index along y-axis
      !! @param [in] iz local index along z-axis
      !! @retval injection on inlet.
      function valueInjection( ix, iy, iz ) result(injection)
         integer, intent(in) :: ix, iy, iz
         integer :: iw
         real(8) :: injection

         iw=lw(ix,iy,iz)
         injection = well(iw)%injection
      end function valueInjection
   end subroutine solvePressure

   !>
   !! @brief calculate velocity.
   subroutine calcVelocity
      !--- Local variables.
      integer :: ix, iy, iz

      !-- Update velocity between grid points.
      do iz=0, param%nz
         do iy=0, param%ny
            do ix=0, param%nx
               ! vx at x-mid and vy at y-mid and vz at z-mid.
               velocity(1:3,ix,iy,iz) = -permeab(1:3,ix,iy,iz)/param%viscosity * grad_pressure(ix,iy,iz)
            end do
         end do
      end do

   contains
      !>
      !! @brief calculate gradient of pressure at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis.
      !! @param [in] iz index along z-axis.
      !! @retval gradient of pressure.
      function grad_pressure( ix, iy, iz ) result(grad)
         integer, intent(in) :: ix, iy, iz
         real(8) :: grad(3)

         grad(1) = (pressure(ix+1,iy,iz)-pressure(ix,iy,iz))/param%dx ! vx at x-mid.
         grad(2) = (pressure(ix,iy+1,iz)-pressure(ix,iy,iz))/param%dx ! vy at y-mid.
         grad(3) = (pressure(ix,iy,iz+1)-pressure(ix,iy,iz))/param%dx ! vy at z-mid.
      end function grad_pressure
   end subroutine calcVelocity

   !>
   !! @brief calculate time derivative of phase.
   subroutine calcPhaseDeriv( taumin )
      use, intrinsic :: ieee_arithmetic  
      !--- Dummy arguments.
      real(8), intent(out) :: taumin
      !--- Local variables.
      integer :: ix, iy, iz
      real(8) :: gpf(3), dpf, normal(3)
      real(8) :: velocity_x, velocity_y, velocity_z, velocity_s, velocity_t
      real(8) :: kxx, tau
      real(8) :: f1, f2, f3, f4, f5, f6, f7
      real(8) :: taumin_mpi

      real(8) :: phix, phiy, phiz, phixx, phiyy, phizz, phixy, phiyz, phizx
      real(8) :: dthetadx, dthetady, dthetadz, dphi_dx, dphi_dy, dphi_dz
      real(8) :: dthetaddpdx, dthetaddpdy, dthetaddpdz, dphi_ddpdx, dphi_ddpdy, dphi_ddpdz
      real(8) :: dwdx, dwdy, dwdz, ddwdphi_x, ddwdphi_y, ddwdphi_z, ddwdthetax, ddwdthetay, ddwdthetaz
      real(8) :: dwdtheta, dwdphi_, ddwdphi_, ddwdtheta, ddwdthetadphi_
      real(8) :: w, w0, norm, theta, theta0, dtheta, phi_, phi_0, dphi_

      real(8) :: eps, c, s, c2, s2, c4, s4, c4phi, s4phi
      real(8) :: maxdphi

      phi_0 = 0
      theta0 = 0
      w0 = 2.0d0 * param%dx

      taumin = 1.0d6 ! set a large value.

      maxdphi = 2.0d0 / (2.0d0 * param%dx)

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx

               gpf(:) = grad_phase(ix,iy,iz)
               dpf    = sqrt(sum(gpf(:)**2)) ! |grad(phase)|

               if( dpf>1.0d-6 ) then
                  normal(:)  = gpf(:)/dpf ! grad(phase)/|grad(phase)|
                  velocity_x = 0.5d0*(velocity(1,ix,iy,iz) + velocity(1,ix-1,iy,iz)) ! vx at grid point.
                  velocity_y = 0.5d0*(velocity(2,ix,iy,iz) + velocity(2,ix,iy-1,iz)) ! vy at grid point.
                  velocity_z = 0.5d0*(velocity(3,ix,iy,iz) + velocity(3,ix,iy,iz-1)) ! vz at grid point.
                  velocity_s = + velocity_x*normal(1) + velocity_y*normal(2) +  velocity_z*normal(3) ! surface component.
                  velocity_t = 0.0d0
               else
                  velocity_s = 0.0d0
                  velocity_t = 0.0d0
               end if

               if( velocity_s>=0.0d0 ) then
                  kxx = param%diffusion
               else
                  kxx = param%diffusion + abs(concent(ix,iy,iz))*param%epsilon &
                     * ( param%ral*velocity_s**2 + param%rat*velocity_t**2 ) &
                     / sqrt(velocity_s**2 + velocity_t**2)
               end if

               tau    = 5.0d0/3.0d0*param%alpha*param%lambda*param%epsilon**2/kxx
               taumin = min(tau,taumin) ! find the minimum of tau.

               if( lw(ix,iy,iz)>0 ) then
                  ! no time derivative at fixed grid point.
                  dphasedt(ix,iy,iz) = 0.0d0
                  cycle
               end if

                              ! --- 一階偏微分 ---
               phix = gpf(1)
               phiy = gpf(2)
               phiz = gpf(3)
               
                              ! --- 二階偏微分 ---
               phixx = (phase(ix+1, iy, iz) - 2.0d0 * phase(ix, iy, iz) + phase(ix-1, iy, iz)) / (param%dx**2)
               phiyy = (phase(ix, iy+1, iz) - 2.0d0 * phase(ix, iy, iz) + phase(ix, iy-1, iz)) / (param%dx**2)
               phizz = (phase(ix, iy, iz+1) - 2.0d0 * phase(ix, iy, iz) + phase(ix, iy, iz-1)) / (param%dx**2)

               phixy = ( phase(ix+1, iy+1, iz) - phase(ix+1, iy-1, iz) &
                     - phase(ix-1, iy+1, iz) + phase(ix-1, iy-1, iz) ) / (4.0d0 * param%dx**2)
               ! ∂^2φ/∂y∂z
               phiyz = (  phase(ix, iy+1, iz+1) - phase(ix, iy-1, iz+1)  &
                     - phase(ix, iy+1, iz-1) + phase(ix, iy-1, iz-1) ) / (4.0d0 * param%dx**2)

               ! ∂^2φ/∂z∂x  (= ∂^2φ/∂x∂z)
               phizx = (  phase(ix+1, iy, iz+1) - phase(ix-1, iy, iz+1)  &
                     - phase(ix+1, iy, iz-1) + phase(ix-1, iy, iz-1) ) / (4.0d0 * param%dx**2)

               if (phix == 0 .and. phiy == 0 .and. phiz == 0 .and. phixx == 0 .and. phiyy == 0 .and. phizz == 0 .and. phixy == 0 .and. phiyz == 0 .and. phizx == 0) then
                  dphasedt(ix,iy,iz) = 0
                  cycle
               endif
               

               theta = atan(phiy / (phix + 1.0d-20))
               phi_ = acos(phiz / (dpf + 1.0d-20))

               dtheta = theta - theta0
               dphi_ = phi_ - phi_0



               ! 前計算（任意：可読性＆計算削減）
               eps   = param%epsilon
               c     = cos(dtheta)
               s     = sin(dtheta)
               c2    = c*c
               s2    = s*s
               c4    = c2*c2
               s4    = s2*s2
               c4phi = cos(4.0d0 * dphi_)
               s4phi = sin(4.0d0 * dphi_)

               ! W(n) = W0( 1 + ε [ 7 cos^4Δθ − 6 cos^2Δθ + sin^4Δθ cos(4Δφ) ] )
               w = w0 * ( 1.0d0 + eps * ( 7.0d0*c4 - 6.0d0*c2 + s4*c4phi ) )

               ! ∂W/∂θ = W0 { 4 ε sinΔθ cosΔθ ( −7 cos^2Δθ + 3 + sin^2Δθ cos(4Δφ) ) }
               dwdtheta = w0 * ( 4.0d0*eps * s * c * ( -7.0d0*c2 + 3.0d0 + s2*c4phi ) )

               ! ∂^2W/∂θ^2 = W0 [ 4 ε { cos(2Δθ)(−7 cos^2Δθ + 3 + sin^2Δθ cos(4Δφ))
               !                         + 2 sin^2Δθ cos^2Δθ (7 + cos(4Δφ)) } ]
               ddwdtheta = w0 * ( 4.0d0*eps * ( cos(2.0d0*dtheta) * ( -7.0d0*c2 + 3.0d0 + s2*c4phi )  &
                              + 2.0d0*s2*c2*( 7.0d0 + c4phi ) ) )

               ! ∂W/∂φ = W0 ( -4 ε sin^4Δθ sin(4Δφ) )
               dwdphi_ = w0 * ( -4.0d0*eps * s4 * s4phi )

               ! ∂^2W/∂φ^2 = W0 ( -16 ε sin^4Δθ cos(4Δφ) )
               ddwdphi_ = w0 * ( -16.0d0*eps * s4 * c4phi )

               ! ∂^2W/∂φ∂θ = W0 ( -16 ε sin^3Δθ cosΔθ sin(4Δφ) )
               ddwdthetadphi_ = w0 * ( -16.0d0*eps * (s2*s) * c * s4phi )
               !  = w0 * ( -16.0d0*eps * sin(dtheta)**3 * cos(dtheta) * sin(4*dphi_) )


               dthetadx = 0.0d0
               dthetady = 0.0d0
               dthetadz = 0.0d0
               if (phix > maxdphi .or. phiy > maxdphi) then
                  dthetadx = ( phixy*phix - phiy*phixx ) / (phix**2 + phiy**2 + 1.0d-20)
                  dthetady = ( phiyy*phix - phiy*phixy ) / (phix**2 + phiy**2 + 1.0d-20)
                  dthetadz = ( phiyz*phix - phiy*phizx ) / (phix**2 + phiy**2 + 1.0d-20)
               endif
               
               dphi_dx = 0
               dphi_dy = 0
               dphi_dz = 0
               if (phix > maxdphi .or. phiy > maxdphi .or. phiz > maxdphi) then
                  dphi_dx = ( phiz * (phix * phixx + phiy * phixy) - (phix**2 + phiy**2) * phizx ) / &
                           ( (phix**2 + phiy**2 + phiz**2) * sqrt(phix**2 + phiy**2) + 1.0d-20)

                  dphi_dy = ( phiz * (phix * phixy + phiy * phiyy) - (phix**2 + phiy**2) * phiyz ) / &
                           ( (phix**2 + phiy**2 + phiz**2) * sqrt(phix**2 + phiy**2) + 1.0d-20)

                  dphi_dz = ( phiz * (phix * phizx + phiy * phiyz) - (phix**2 + phiy**2) * phizz ) / &
                           ( (phix**2 + phiy**2 + phiz**2) * sqrt(phix**2 + phiy**2) + 1.0d-20)
               endif

               dthetaddpdx = 0.0d0
               dthetaddpdy = 0.0d0
               if (phix > maxdphi .or. phiy > maxdphi) then
                  dthetaddpdx = - phiy * (phix**2 + phiy**2 + phiz**2) / ((phix**2 + phiy**2) + 1.0d-20)
                  dthetaddpdy =   phix * (phix**2 + phiy**2 + phiz**2) / ((phix**2 + phiy**2) + 1.0d-20)
               endif
               dthetaddpdz = 0.0d0

               ! --- 各方向の ∂φ/∂(∂φ/∂x), ∂φ/∂(∂φ/∂y), ∂φ/∂(∂φ/∂z) の導関数 ---

               dphi_ddpdx =   phix * phiz / (sqrt(phix**2 + phiy**2) * dpf + 1.0d-20)
               dphi_ddpdy =   phiy * phiz / (sqrt(phix**2 + phiy**2) * dpf + 1.0d-20)
               dphi_ddpdz = - sqrt(phix**2 + phiy**2) / (dpf + 1.0d-20)

               dwdx = dwdtheta * dthetadx + dwdphi_ * dphi_dx
               dwdy = dwdtheta * dthetady + dwdphi_ * dphi_dy
               dwdz = dwdtheta * dthetadz + dwdphi_ * dphi_dz

               ddwdthetax = ddwdtheta * dthetadx + ddwdthetadphi_ * dphi_dx
               ddwdthetay = ddwdtheta * dthetady + ddwdthetadphi_ * dphi_dy
               ddwdthetaz = ddwdtheta * dthetadz + ddwdthetadphi_ * dphi_dz

               ddwdphi_x = ddwdphi_ * dphi_dx + ddwdthetadphi_ * dthetadx
               ddwdphi_y = ddwdphi_ * dphi_dy + ddwdthetadphi_ * dthetady
               ddwdphi_z = ddwdphi_ * dphi_dz + ddwdthetadphi_ * dthetadz      

               

               f1    = (1.0d0-phase(ix,iy,iz)**2) &
                  * (phase(ix,iy,iz)-param%lambda*concent(ix,iy,iz))

               
               ! ---- 1) 2W (∇W · ∇φ)
               f2 = 2.0d0 * w * ( dwdx * phix + dwdy * phiy + dwdz * phiz )
               ! f2 = 0

               ! ---- 2) W^2 ∇^2φ
               f3 = w**2 * lap_phase(ix, iy, iz)

               ! ---- 3) 角依存ブロック（x偏微分の鎖則部分）
               !   ∂W/∂x { (∂W/∂θ) * ∂θ/∂(∂φ/∂x) * |∇φ|^2 + (∂W/∂ϕ) * ∂ϕ/∂(∂φ/∂x) * |∇φ|^2 }
               ! + W { ∂^2W/∂θ∂x * ∂θ/∂(∂φ/∂x) * |∇φ|^2 + ∂^2W/∂ϕ∂x * ∂ϕ/∂(∂φ/∂x) * |∇φ|^2 }
               f4 = dwdx * ( dwdtheta * dthetaddpdx * (dpf**2) + dwdphi_ * dphi_ddpdx * (dpf**2) )  &
                  + w    * ( ddwdthetax * dthetaddpdx * (dpf**2) + ddwdphi_x * dphi_ddpdx * (dpf**2) )
               ! f4 = 0

               ! ---- 4) 角依存ブロック（y偏微分の鎖則部分）
               !   ∂W/∂y {... x と同様 ...}
               f5 = dwdy * ( dwdtheta * dthetaddpdy * (dpf**2) + dwdphi_ * dphi_ddpdy * (dpf**2) )  &
                  + w    * ( ddwdthetay * dthetaddpdy * (dpf**2) + ddwdphi_y * dphi_ddpdy * (dpf**2) )
               ! f5 = 0

               ! ---- 5) 角依存ブロック（z偏微分の鎖則部分）
               !   ∂W/∂z {... x と同様 ...}
               f6 = dwdz * ( dwdtheta * dthetaddpdz * (dpf**2) + dwdphi_ * dphi_ddpdz * (dpf**2) )  &
                  + w    * ( ddwdthetaz * dthetaddpdz * (dpf**2) + ddwdphi_z * dphi_ddpdz * (dpf**2) )
               ! f6 = 0

               !f7    = -param%epsilon**2* div_normal(ix,iy,iz) * dpf
               f7 = 0

               dphasedt(ix,iy,iz) = (f1+f2+f3+f4+f5+f6+f7)/tau

                  ! if (velocity_s<0.0d0) then
                  !    write(*,*)
                  ! endif


               ! ! if (ieee_is_nan(dphasedt(ix,iy,iz))) then
               ! !    write(*,*)
               ! ! endif

                  if (phase(ix,iy,iz) > 1 .or. phase(ix,iy,iz) < -1) then
                     write(*,*)"aa"
                  endif
               

               ! if (ix == 5 .and. iy == 4 .and. iz == 4) then
               !    write(*,*)
               ! endif


            end do
         end do
      end do

      taumin_mpi = taumin
      call MPI_Allreduce( taumin_mpi, taumin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_ierr )

   contains


      !>
      !! @brief calculate laplacian of phase at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis
      !! @param [in] iz index along z-axis
      !! @retval laplacian of phase.
      function lap_phase( ix, iy, iz ) result(lap)
         integer, intent(in) :: ix, iy, iz
         real(8) :: lap

         lap = ( &
            +    phase(ix+1,iy,iz) + phase(ix-1,iy,iz) &
            +    phase(ix,iy+1,iz) + phase(ix,iy-1,iz) &
            +    phase(ix,iy,iz+1) + phase(ix,iy,iz-1) &
            -  6*phase(ix,iy,iz) )/param%dx**2
      end function lap_phase

      !>
      !! @brief calculate gradient of phase at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis
      !! @param [in] iz index along z-axis
      !! @retval gradient of phase.
      function grad_phase( ix, iy, iz ) result(grad)
         integer, intent(in) :: ix, iy, iz
         real(8) :: grad(3)

         grad(1) = ( phase(ix+1,iy,iz) - phase(ix-1,iy,iz) )/(2.0d0*param%dx)
         grad(2) = ( phase(ix,iy+1,iz) - phase(ix,iy-1,iz) )/(2.0d0*param%dx)
         grad(3) = ( phase(ix,iy,iz+1) - phase(ix,iy,iz-1) )/(2.0d0*param%dx)
      end function grad_phase

      !>
      !! @brief calculate divergence of normal vector at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis
      !! @param [in] iz index along z-axis
      !! @retval divergence of normal vector.
      function div_normal( ix, iy, iz ) result(div)
         integer, intent(in) :: ix, iy, iz
         real(8) :: div

         !--- Local variables.
         real(8) :: phase_ppc_cpc, phase_pmc_cmc, phase_mpc_cpc, phase_mmc_cmc ! phases at middle points
         real(8) :: phase_ppc_pcc, phase_mpc_mcc, phase_pmc_pcc, phase_mmc_mcc ! phases at middle points
         real(8) :: phase_pcp_ccp, phase_pcm_ccm, phase_mcp_ccp, phase_mcm_ccm ! phases at middle points
         real(8) :: phase_cpp_ccp, phase_cpm_ccm, phase_cmp_ccp, phase_cmm_ccm ! phases at middle points
         real(8) :: phase_pcp_pcc, phase_mcp_mcc, phase_cpp_cpc, phase_cmp_cmc ! phases at middle points
         real(8) :: phase_pcc_pcm, phase_mcc_mcm, phase_cpc_cpm, phase_cmc_cmm ! phases at middle points
         real(8) :: gphase_xp(3), gphase_xm(3), gphase_yp(3), gphase_ym(3), gphase_zp(3), gphase_zm(3) ! gradient of phase 
         real(8) :: normal_xp, normal_xm, normal_yp, normal_ym, normal_zp, normal_zm

         phase_ppc_cpc = ( phase(ix+1,iy+1,iz) + phase(ix,iy+1,iz) )*0.5d0 ! phase(ix+1/2, iy+1, iz)
         phase_pmc_cmc = ( phase(ix+1,iy-1,iz) + phase(ix,iy-1,iz) )*0.5d0 ! phase(ix+1/2, iy-1, iz)
         phase_mpc_cpc = ( phase(ix-1,iy+1,iz) + phase(ix,iy+1,iz) )*0.5d0 ! phase(ix-1/2, iy+1, iz)
         phase_mmc_cmc = ( phase(ix-1,iy-1,iz) + phase(ix,iy-1,iz) )*0.5d0 ! phase(ix-1/2, iy-1, iz)
         phase_ppc_pcc = ( phase(ix+1,iy+1,iz) + phase(ix+1,iy,iz) )*0.5d0 ! phase(ix+1, iy+1/2, iz)
         phase_mpc_mcc = ( phase(ix-1,iy+1,iz) + phase(ix-1,iy,iz) )*0.5d0 ! phase(ix-1, iy+1/2, iz)
         phase_pmc_pcc = ( phase(ix+1,iy-1,iz) + phase(ix+1,iy,iz) )*0.5d0 ! phase(ix+1, iy-1/2, iz)
         phase_mmc_mcc = ( phase(ix-1,iy-1,iz) + phase(ix-1,iy,iz) )*0.5d0 ! phase(ix-1, iy-1/2, iz)
         phase_pcp_ccp = ( phase(ix+1,iy,iz+1) + phase(ix,iy,iz+1) )*0.5d0 ! phase(ix+1/2, iy, iz+1)
         phase_pcm_ccm = ( phase(ix+1,iy,iz-1) + phase(ix,iy,iz-1) )*0.5d0 ! phase(ix+1/2, iy, iz-1)
         phase_mcp_ccp = ( phase(ix-1,iy,iz+1) + phase(ix,iy,iz+1) )*0.5d0 ! phase(ix-1/2, iy, iz+1)
         phase_mcm_ccm = ( phase(ix-1,iy,iz-1) + phase(ix,iy,iz-1) )*0.5d0 ! phase(ix-1/2, iy, iz-1)
         phase_cpp_ccp = ( phase(ix,iy+1,iz+1) + phase(ix,iy,iz+1) )*0.5d0 ! phase(ix, iy+1/2, iz+1)
         phase_cpm_ccm = ( phase(ix,iy+1,iz-1) + phase(ix,iy,iz-1) )*0.5d0 ! phase(ix, iy+1/2, iz-1)
         phase_cmp_ccp = ( phase(ix,iy-1,iz+1) + phase(ix,iy,iz+1) )*0.5d0 ! phase(ix, iy-1/2, iz+1)
         phase_cmm_ccm = ( phase(ix,iy-1,iz-1) + phase(ix,iy,iz-1) )*0.5d0 ! phase(ix, iy-1/2, iz-1)
         phase_pcp_pcc = ( phase(ix+1,iy,iz+1) + phase(ix+1,iy,iz) )*0.5d0 ! phase(ix+1, iy, iz+1/2)
         phase_mcp_mcc = ( phase(ix-1,iy,iz+1) + phase(ix-1,iy,iz) )*0.5d0 ! phase(ix-1, iy, iz+1/2)
         phase_cpp_cpc = ( phase(ix,iy+1,iz+1) + phase(ix,iy+1,iz) )*0.5d0 ! phase(ix, iy+1, iz+1/2)
         phase_cmp_cmc = ( phase(ix,iy-1,iz+1) + phase(ix,iy-1,iz) )*0.5d0 ! phase(ix, iy-1, iz+1/2)
         phase_pcc_pcm = ( phase(ix+1,iy,iz) + phase(ix+1,iy,iz-1) )*0.5d0 ! phase(ix+1, iy, iz-1/2)
         phase_mcc_mcm = ( phase(ix-1,iy,iz) + phase(ix-1,iy,iz-1) )*0.5d0 ! phase(ix-1, iy, iz-1/2)
         phase_cpc_cpm = ( phase(ix,iy+1,iz) + phase(ix,iy+1,iz-1) )*0.5d0 ! phase(ix, iy+1, iz-1/2)
         phase_cmc_cmm = ( phase(ix,iy-1,iz) + phase(ix,iy-1,iz-1) )*0.5d0 ! phase(ix, iy-1, iz-1/2)

         gphase_xp(1) = ( phase(ix+1,iy,iz) - phase(ix,iy,iz))    /param%dx     ! grad_x phase(ix+1/2,iy,iz)
         gphase_xp(2) = ( phase_ppc_cpc     - phase_pmc_cmc )     /(2*param%dx) ! grad_y phase(ix+1/2,iy,iz)
         gphase_xp(3) = ( phase_pcp_ccp     - phase_pcm_ccm )     /(2*param%dx) ! grad_z phase(ix+1/2,iy,iz)
         gphase_xm(1) = ( phase(ix,iy,iz)   - phase(ix-1,iy,iz) ) /param%dx     ! grad_x phase(ix-1/2,iy,iz)
         gphase_xm(2) = ( phase_mpc_cpc     - phase_mmc_cmc )     /(2*param%dx) ! grad_y phase(ix-1/2,iy,iz)
         gphase_xm(3) = ( phase_mcp_ccp     - phase_mcm_ccm )     /(2*param%dx) ! grad_z phase(ix-1/2,iy,iz)
         gphase_yp(1) = ( phase_ppc_pcc     - phase_mpc_mcc )     /(2*param%dx) ! grad_x phase(ix,iy+1/2,iz)
         gphase_yp(2) = ( phase(ix,iy+1,iz) - phase(ix,iy,iz) )   /param%dx     ! grad_y phase(ix,iy+1/2,iz)
         gphase_yp(3) = ( phase_cpp_ccp     - phase_cpm_ccm )     /(2*param%dx) ! grad_z phase(ix,iy+1/2,iz)
         gphase_ym(1) = ( phase_pmc_pcc     - phase_mmc_mcc )     /(2*param%dx) ! grad_x phase(ix,iy-1/2,iz)
         gphase_ym(2) = ( phase(ix,iy,iz)   - phase(ix,iy-1,iz) ) /param%dx     ! grad_y phase(ix,iy-1/2,iz)
         gphase_ym(3) = ( phase_cmp_ccp     - phase_cmm_ccm )     /(2*param%dx) ! grad_z phase(ix,iy-1/2,iz)
         gphase_zp(1) = ( phase_pcp_pcc     - phase_mcp_mcc )     /(2*param%dx) ! grad_x phase(ix,iy,iz+1/2)
         gphase_zp(2) = ( phase_cpp_cpc     - phase_cmp_cmc )     /(2*param%dx) ! grad_y phase(ix,iy,iz+1/2)
         gphase_zp(3) = ( phase(ix,iy,iz+1) - phase(ix,iy,iz) )   /param%dx     ! grad_z phase(ix,iy,iz+1/2)
         gphase_zm(1) = ( phase_pcc_pcm     - phase_mcc_mcm )     /(2*param%dx) ! grad_x phase(ix,iy,iz-1/2)
         gphase_zm(2) = ( phase_cpc_cpm     - phase_cmc_cmm )     /(2*param%dx) ! grad_y phase(ix,iy,iz-1/2)
         gphase_zm(3) = ( phase(ix,iy,iz)   - phase(ix,iy,iz-1) ) /param%dx     ! grad_z phase(ix,iy,iz-1/2)

         normal_xp  = merge( 0.0d0, gphase_xp(1)/sqrt(sum(gphase_xp(:)**2)), gphase_xp(1)==0.0d0 ) ! n_x(ix+1/2,iy,iz)
         normal_xm  = merge( 0.0d0, gphase_xm(1)/sqrt(sum(gphase_xm(:)**2)), gphase_xm(1)==0.0d0 ) ! n_x(ix-1/2,iy,iz)
         normal_yp  = merge( 0.0d0, gphase_yp(2)/sqrt(sum(gphase_yp(:)**2)), gphase_yp(2)==0.0d0 ) ! n_y(ix,iy+1/2,iz)
         normal_ym  = merge( 0.0d0, gphase_ym(2)/sqrt(sum(gphase_ym(:)**2)), gphase_ym(2)==0.0d0 ) ! n_y(ix,iy-1/2,iz)
         normal_zp  = merge( 0.0d0, gphase_zp(3)/sqrt(sum(gphase_zp(:)**2)), gphase_zp(3)==0.0d0 ) ! n_z(ix,iy,iz+1/2)
         normal_zm  = merge( 0.0d0, gphase_zm(3)/sqrt(sum(gphase_zm(:)**2)), gphase_zm(3)==0.0d0 ) ! n_z(ix,iy,iz-1/2)

         div = (normal_xp-normal_xm)/param%dx + (normal_yp-normal_ym)/param%dx + (normal_zp-normal_zm)/param%dx
      end function div_normal
   end subroutine calcPhaseDeriv

   !>
   !! @brief calculate time derivative of concentration.
   subroutine calcConcentDeriv
      !--- Local variables.
      integer :: ix, iy, iz
      real(8) :: c1, c2, c3, c7

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               if( lw(ix,iy,iz)>0 ) then
                  ! no time derivative at fixed grid point.
                  dconcentdt(ix,iy,iz) = 0.0d0
                  cycle
               end if

               c1 = param%diffusion * lap_concent(ix,iy,iz)
               c2 = param%alpha * dphasedt(ix,iy,iz)
               c3 = 0.0d0
               c7 = - div_flow(ix,iy,iz)

               dconcentdt(ix,iy,iz) = c1+c2+c3+c7
            end do
         end do
      end do

   contains

      !>
      !! @brief calculate divergence of flow at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis
      !! @param [in] iz index along z-axis
      !! @retval divergence of flow.
      function div_flow( ix, iy, iz ) result(div)
         integer, intent(in) :: ix, iy, iz
         real(8) :: div
         real(8) :: flow_xp, flow_xm, flow_yp, flow_ym, flow_zp, flow_zm

         ! flow as velocity * concent using upwind.
         flow_xp = velocity(1,ix,iy,iz) &
            * merge( concent(ix,iy,iz),   concent(ix+1,iy,iz), velocity(1,ix,iy,iz)  >0.0d0 )
         flow_xm = velocity(1,ix-1,iy,iz) &
            * merge( concent(ix-1,iy,iz), concent(ix,iy,iz),   velocity(1,ix-1,iy,iz)>0.0d0 )
         flow_yp = velocity(2,ix,iy,iz) &
            * merge( concent(ix,iy,iz),   concent(ix,iy+1,iz), velocity(2,ix,iy,iz)  >0.0d0 )
         flow_ym = velocity(2,ix,iy-1,iz) &
            * merge( concent(ix,iy-1,iz), concent(ix,iy,iz),   velocity(2,ix,iy-1,iz)>0.0d0 )
         flow_zp = velocity(3,ix,iy,iz) &
            * merge( concent(ix,iy,iz),   concent(ix,iy,iz+1), velocity(3,ix,iy,iz)  >0.0d0 )
         flow_zm = velocity(3,ix,iy,iz-1) &
            * merge( concent(ix,iy,iz-1), concent(ix,iy,iz),   velocity(3,ix,iy,iz-1)>0.0d0 )

         div = (flow_xp-flow_xm)/param%dx + (flow_yp-flow_ym)/param%dx + (flow_zp-flow_zm)/param%dx
      end function div_flow

      !>
      !! @brief calculate laplacian of concentration at a given grid point.
      !! @param [in] ix index along x-axis.
      !! @param [in] iy index along y-axis
      !! @param [in] iz index along z-axis      
      !! @retval laplacian of concentration.
      function lap_concent( ix, iy, iz ) result(lap)
         integer, intent(in) :: ix, iy, iz
         real(8) :: lap

         lap = ( &
            +    concent(ix+1,iy,iz) + concent(ix-1,iy,iz) &
            +    concent(ix,iy+1,iz) + concent(ix,iy-1,iz) &
            +    concent(ix,iy,iz+1) + concent(ix,iy,iz-1) &
            -  6*concent(ix,iy,iz) )/param%dx**2
      end function lap_concent
   end subroutine calcConcentDeriv

   !>
   !! @brief evolve phase.
   subroutine evolvePhase
      phase(:,:,:) = phase(:,:,:) + dphasedt(:,:,:)*param%dtime

      call updateHalo( phase )
   end subroutine evolvePhase

   !>
   !! @brief evolve concentration.
   subroutine evolveConcent
      concent(:,:,:) = concent(:,:,:) + dconcentdt(:,:,:)*param%dtime

      call updateHalo( concent )
   end subroutine evolveConcent

   !>
   !! @brief update field at halo grid points.
   !! @param field [in,out] field whose values at halo grid points to be update.
   subroutine updateHalo( field )
      !--- Dummy arguments.
      real(8), intent(inout) :: field(0:param%nx+1,0:param%ny+1,0:param%nz+1)

      call MPI__SendrecvHalo( field )

      if( param%ex(1) ) field(0,:,:)          = field(1,:,:)          ! x- halo
      if( param%ex(2) ) field(param%nx+1,:,:) = field(param%nx+0,:,:) ! x+ halo
      if( param%ey(1) ) field(:,0,:)          = field(:,1,:)          ! y- halo
      if( param%ey(2) ) field(:,param%ny+1,:) = field(:,param%ny+0,:) ! y+ halo
      if( param%ez(1) ) field(:,:,0)          = field(:,:,1)          ! z- halo
      if( param%ez(2) ) field(:,:,param%nz+1) = field(:,:,param%nz+0) ! z+ halo
   end subroutine updateHalo

   !>
   !! @brief show progress to screen and data file.
   !! @param [in] filename file name of output data file.
   !! @param [in] time current time.
   !! @param [in] itime current time step.
   !! @param [in] pinj  injection pressure on inlet.
   subroutine showProgress( filename, time, itime, pinj )
      !--- Dummy arguments.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: itime
      real(8), intent(in) :: time
      real(8), intent(in) :: pinj
      !--- Local variables.
      integer :: ix, iy, iz
      integer :: iw
      real(8) :: x, y, z, xref, yref, zref, r2, r2whmax, rwhmax
      real(8) :: vcal, vacid, vcaco3, qinj, b, kc, xchl
      real(8) :: phase_max, phase_min, concent_max, concent_min
      real(8) :: phase_max_mpi, phase_min_mpi, concent_max_mpi, concent_min_mpi
      real(8) :: pv_mpi, vcal_mpi, rwhmax_mpi
      logical, save :: firsttime = .true.
      integer, save :: ounit
      real(8), save :: pv
      real(8), save :: cpu_tstr, cpu_tend, cpu_ttot
      integer(8), save :: time_begin_c, time_end_c, countpersec, countmax

      xchl = param%cwtp*param%b100*param%rhoa/param%rhos

      if( firsttime ) then
         if( mpi_parent ) open( newunit=ounit, file=filename, status="replace" )

         !--- Record start times.
         call system_clock( time_begin_c, countpersec, countMax )
         call CPU_TIME( cpu_tstr )

         !--- Compute initial pore volume of rock.
         pv = param%dx**3 * param%porosity * 0.5d0*sum(1.0d0-phase(1:param%nx,1:param%ny,1:param%nz))
         pv_mpi = pv
         call MPI_Allreduce( pv_mpi, pv,1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr )

         b  = param%rhoa/param%rhos
         kc = param%cwtp*param%b100/(1.0d0-param%porosity)

         !--- write input data in output file.
         if( mpi_parent ) then
            write(ounit,*) "***** CALCULATION SUMMARY FOR PFACID PROGRAM *****"
            write(ounit,*) "/// Input data ///"
            write(ounit,*) "*** Grid data"
            write(ounit,"(3x,'Num of grids, nx,ny,nz =',3i6)") param%nx_g, param%ny_g, param%nz_g
            write(ounit,"(3x,'Grid size, dx       =',es14.6,' [m]')") param%dx
            write(ounit,"(3x,'Core length, Lx, Ly, Lz =',3f10.4,' [m]')") param%lx, param%ly, param%lz
            write(ounit,*) "*** Rock and fluid data"
            write(ounit,"(3x,'acid con., cwtp      =',f12.4,' [kg HCl/kg solution]')") param%cwtp
            write(ounit,"(3x,'acid density, rhoa   =',f12.1,' [kg/m3]')") param%rhoa
            write(ounit,"(3x,'rock density, rhos   =',f12.1,' [kg/m3]')") param%rhos
            write(ounit,"(3x,'viscosity, visco     =',es14.6,' [Pa*s]')") param%viscosity
            write(ounit,"(3x,'rock porosity, poro  =',f12.6,' [dec]')") param%porosity
            write(ounit,"(3x,'diffusivity coef. ,d =',es12.4,' [m2/s]')") param%diffusion
            write(ounit,"(3x,'gravimetric dissolving power, b100 =',f8.3,' [kg CaCO3/kg HCl]')") param%b100
            write(ounit,"(3x,'volumetric dissolving power, Xacid =',f8.3,' [m3 CaCO3/m3 HCl]')") xchl
            write(ounit,*) "*** Phase-field related parameters"
            write(ounit,"(3x,'lambda =',f12.6,' [-]')") param%lambda
            write(ounit,"(3x,'kc     =',f12.6,' [-]')") kc
            write(ounit,"(3x,'b      =',f12.6,' [-]')") b
            write(ounit,"(3x,'alpha  =',f12.6,' [-]')") param%alpha
            write(ounit,*) "*** Injection data"
            write(ounit,"(3x,'time step size, dt     =',es12.4,' [s]')") param%dtime
            write(ounit,"(3x,'initial pressure, pini =',es12.4,' [m/s]')") param%pini

            write(*,"('*** 3D Phase-Field Acid Wormholing Model, wPFAW3D ***')")
            write(*,"(/1x,'# of equations,     neq =',i12)") param%ngrid_g
            ! write(*,"(1x,'# of nonzero entries, m =',i8)") nnz     ! Symmetric
            write(*,"(1x,'pore volume =',es20.6,' [m3]')") pv
            write(ounit,"(1x,'pore volume =',es20.6,' [m3]')") pv
            write(*,"(1x,'Initial time step size =',es12.4,' [s]')") param%dtime
            write(*,"(1x,'Start time =',es12.4,' [s]')") param%time_start
            write(*,"(1x,'End time   =',es12.4,' [s]')") param%time_end
            !write(*,"(1x,'Total Injection time =',es12.4,' [s]')") dt*ntime
            write(*,*)
            write(*,"(1x,'  ITIME     TIME        PVINJ     MAX_PHI MIN_PHI  MAX_C   MIN_C     PINJ         RWH         DT   ')")
            write(*,"(1x,'   (-)       (s)        (v/v)       (-)     (-)     (-)     (-)      (Pa)         (m)         (s)  ')")
            write(*,"(1x,' ------  -----------  ---------   ------- ------- ------- ------- -----------  ----------  --------')")

            !--- header lines for output file.
            write(ounit,"(//'/// SIMULATION RESULTS ///')")
            write(ounit,"(13x,'ITIME',2x,'TIME',15x,'PVINJ',13x,'WELL_PRESS'," // &
                 "10x,'WH_LENGTH',10x,'INJ_ACID',11x,'DISSOLVED_CaCO3',10x,'DISS.CaCO3_c')")
            write(ounit,"(13x,'(-)',4x,'(sec)',15x,'(v/v)',15x,'(Pa)'," // &
                 "16x,'(m)',14x,'(m3/m)',14x,'(m3/m)',14x,'(m3/m)')")
            write(ounit,"(5x,'-----',2x,6('   --------------   '))")
         end if

         firsttime = .false.
      end if

      !--- Simulation results.
      ! compute max wormhole length.
      r2whmax = 0.0d0
      xref = param%dx*(param%ix_ref-1) + 0.5d0*param%dx
      yref = param%dx*(param%iy_ref-1) + 0.5d0*param%dx
      zref = param%dx*(param%iz_ref-1) + 0.5d0*param%dx

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               if( phase(ix,iy,iz)>0.0d0 ) then
                  x  = param%dx*(ix+param%sx-1) + 0.5d0*param%dx
                  y  = param%dx*(iy+param%sy-1) + 0.5d0*param%dx
                  z  = param%dx*(iz+param%sz-1) + 0.5d0*param%dx
                  r2 = (x-xref)**2+(y-yref)**2+(z-zref)**2
                  r2whmax = max(r2whmax,r2)
               end if
            end do
         end do
      end do

      rwhmax = sqrt(r2whmax)
      vcal = param%dx**3 * (1.0d0-param%porosity) * 0.5d0*sum(1.0d0+phase(1:param%nx,1:param%ny,1:param%nz))

      rwhmax_mpi = rwhmax
      call MPI_Allreduce( rwhmax_mpi, rwhmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr )
      vcal_mpi = vcal
      call MPI_Allreduce( vcal_mpi, vcal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr )

      !--- Compute total injection rate in m3/m.
      qinj=0.0d0
      do iw=1, param%nwell
         if( well(iw)%injection/=0.0d0 ) qinj=qinj+well(iw)%injection
      end do

      vacid  = qinj*time
      vcaco3 = qinj*time*xchl

      if( mpi_parent ) then
         write(ounit, '(1x,i8,1x,7es20.10)') &
            itime, time, vacid/pv, &
            pinj, rwhmax, vacid, vcaco3, vcal
      end if

      phase_max   = maxval(phase)
      phase_min   = minval(phase)
      concent_max = maxval(concent)
      concent_min = minval(concent)

      phase_max_mpi = phase_max
      call MPI_Allreduce( phase_max_mpi, phase_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr )
      phase_min_mpi = phase_min
      call MPI_Allreduce( phase_min_mpi, phase_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_ierr )
      concent_max_mpi = concent_max
      call MPI_Allreduce( concent_max_mpi, concent_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr )
      concent_min_mpi = concent_min
      call MPI_Allreduce( concent_min_mpi, concent_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_ierr )

      if( mpi_parent ) then
         write(*,"(i8,es13.5,es11.3,1x,4(f8.3),es13.5,es12.4,es10.2)") &
            itime, time, vacid/pv, &
            phase_max, phase_min, concent_max, concent_min, pinj, rwhmax, param%dtime
      end if

      return

      entry doneProgress
      
      call system_clock( time_end_c )
      call cpu_time( cpu_tend )
      cpu_ttot = (cpu_tend-cpu_tstr)/60.0d0

      if( mpi_parent ) then
         write(ounit,"(/'CPU time =',es12.4,' [min]')") cpu_ttot
         write(ounit,"('Cal time =',es12.4,' [min]')") dble(time_end_c - time_begin_c)/CountPerSec/60.0d0
         write(*,"(/'CPU time =',es12.4,' [min]')") cpu_ttot
         write(*,"('Cal time =',es12.4,' [min]')") dble(time_end_c - time_begin_c)/CountPerSec/60.0d0

         close(ounit)
      end if

      return
   end subroutine showProgress

   !>
   !! @brief save fields in result file.
   !! @param [in] filename file name of output result file.
   !! @param [in] time current time.
   subroutine saveResult( filename, time )
      !--- Dummy arguments.
      character(len=*), intent(in) :: filename
      real(8), intent(in) :: time
      !--- Local variables.
      logical, save :: firsttime = .true.
      integer, save :: ounit
      integer :: status
      integer(8) :: offset
      integer, parameter :: ByteOrderCheck = int(Z'091d')
      real, parameter :: zeros(3) = [0.0,0.0,0.0]
      integer :: ix_g, iy_g, iz_g
      integer :: ix, iy, iz
      character(len=1024) :: line

      integer :: ngrid_each(0:mpi_size-1)
      integer :: ngrid_before, ngrid_after

      type GiDcoord_type
         integer :: index
         real :: x, y, z
      end type GiDcoord_Type
      type(GiDcoord_type) :: gidcoord(param%nx_g+1)

      type GiDelement_type
         integer :: index
         integer :: nid(8)
         integer :: conn
      end type GiDelement_Type
      type(GiDelement_type) :: gidelement(param%nx_g)

      type GiDvalue_type
         integer :: index
         real :: value
      end type GiDvalue_Type
      type(GiDvalue_type) :: gidvalue(param%nx)

      !--- Output mesh data in GiD format only at the first time.
      if( firsttime ) then
         call MPI_File_delete( filename, MPI_INFO_NULL, mpi_ierr )
         call MPI_File_open( MPI_COMM_WORLD, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, ounit, status )
         if( status/=0 ) then
            write(*,*) "Error: can not create result file: ", trim(filename)
            call exit(1)
         end if

         call MPI_File_write_all( ounit, ByteOrderCheck, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )
         call write_string('GiDPostEx1.1')

         call write_string('MESH "quadmesh" dimension 3 ElemType Hexahedra Nnode 8')
         call write_string('# color 0.7 0.7 0.4')
         call write_string('Unit "m"')
         call write_string('Coordinates -1 Indexed')

         do iz_g=1, param%nz_g+1
            do iy_g=1, param%ny_g+1
               do ix_g=1, param%nx_g+1
                  gidcoord(ix_g)%index = mindex3D(ix_g,iy_g,iz_g)
                  gidcoord(ix_g)%x     = real((ix_g-1)*param%dx,kind=4)
                  gidcoord(ix_g)%y     = real((iy_g-1)*param%dx,kind=4)
                  gidcoord(ix_g)%z     = real((iz_g-1)*param%dx,kind=4)
               end do
               call MPI_File_write_all( ounit, gidcoord, int(sizeof(gidcoord)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
            end do
         end do
         call write_string('End Coordinates')

         call write_string('Elements -1 Indexed')

         do iz_g=1, param%nz_g
            do iy_g=1, param%ny_g
               do ix_g=1, param%nx_g
                  gidelement(ix_g)%index  = eindex3D(ix_g,iy_g,iz_g)
                  gidelement(ix_g)%nid(1) = mindex3D(ix_g+0,iy_g+0,iz_g+0)
                  gidelement(ix_g)%nid(2) = mindex3D(ix_g+1,iy_g+0,iz_g+0)
                  gidelement(ix_g)%nid(3) = mindex3D(ix_g+1,iy_g+1,iz_g+0)
                  gidelement(ix_g)%nid(4) = mindex3D(ix_g+0,iy_g+1,iz_g+0)
                  gidelement(ix_g)%nid(5) = mindex3D(ix_g+0,iy_g+0,iz_g+1)
                  gidelement(ix_g)%nid(6) = mindex3D(ix_g+1,iy_g+0,iz_g+1)
                  gidelement(ix_g)%nid(7) = mindex3D(ix_g+1,iy_g+1,iz_g+1)
                  gidelement(ix_g)%nid(8) = mindex3D(ix_g+0,iy_g+1,iz_g+1)
                  gidelement(ix_g)%conn   = 1
               end do
               call MPI_File_write_all( ounit, gidelement, int(sizeof(gidelement)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
            end do
         end do
         call write_string('End Elements')

         call write_string('GaussPoints "Cell_Center" ElemType Hexahedra')
         call write_string('Number Of Gauss Points: 1')
         call write_string('Natural Coordinates: Given')
         call MPI_File_write_all( ounit, zeros, 3, MPI_FLOAT, MPI_STATUS_IGNORE, mpi_ierr )   
         call write_string('End GaussPoints')

         call MPI_File_sync( ounit, mpi_ierr )

         firsttime = .false.
      end if ! firsttime

      call MPI_Allgather( param%ngrid, 1, MPI_INTEGER, ngrid_each, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr )
      ngrid_before = sum(ngrid_each(0:mpi_rank-1))
      ngrid_after  = sum(ngrid_each(mpi_rank+1:mpi_size-1))

      !--- Output pressure field in GiD format.
      write(line,"('Result ""Pressure"" ""Analysis"" ',es12.5,' Scalar OnGaussPoints ""Cell_Center""')") time
      call write_string(line)
      call write_string('Unit "Pa"')
      call write_string('Values -1 Indexed')

      offset = ngrid_before * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               gidvalue(ix)%index = eindex3D(ix+param%sx,iy+param%sy,iz+param%sz)
               gidvalue(ix)%value = real(pressure(ix,iy,iz),kind=4)
            end do
            call MPI_File_write( ounit, gidvalue, int(sizeof(gidvalue)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
         end do
      end do

      offset = ngrid_after * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      call MPI_File_write_all( ounit, -1, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )   
      call write_string('End Values')

      !--- Output phase field in GiD format.
      write(line,"('Result ""Phase-field"" ""Analysis"" ',es12.5,' Scalar OnGaussPoints ""Cell_Center""')") time
      call write_string(line)
      call write_string('Unit "-"')
      call write_string('Values -1 Indexed')

      offset = ngrid_before * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               gidvalue(ix)%index = eindex3D(ix+param%sx,iy+param%sy,iz+param%sz)
               gidvalue(ix)%value = real(phase(ix,iy,iz),kind=4)
            end do
            call MPI_File_write( ounit, gidvalue, int(sizeof(gidvalue)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
         end do
      end do

      offset = ngrid_after * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      call MPI_File_write_all( ounit, -1, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )   
      call write_string('End Values')

      !--- Output concentration field in GiD format.
      write(line,"('Result ""Concentration"" ""Analysis"" ',es12.5,' Scalar OnGaussPoints ""Cell_Center""')") time
      call write_string(line)
      call write_string('Unit "-"')
      call write_string('Values -1 Indexed')

      offset = ngrid_before * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               gidvalue(ix)%index = eindex3D(ix+param%sx,iy+param%sy,iz+param%sz)
               gidvalue(ix)%value = real(concent(ix,iy,iz),kind=4)
            end do
            call MPI_File_write( ounit, gidvalue, int(sizeof(gidvalue)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
         end do
      end do

      offset = ngrid_after * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      call MPI_File_write_all( ounit, -1, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )   
      call write_string('End Values')

      !--- Output permeability field in GiD format.
      write(line,"('Result ""Permeability"" ""Analysis"" ',es12.5,' Scalar OnGaussPoints ""Cell_Center""')") time
      call write_string(line)
      call write_string('Unit "-"')
      call write_string('Values -1 Indexed')

      offset = ngrid_before * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      do iz=1, param%nz
         do iy=1, param%ny
            do ix=1, param%nx
               gidvalue(ix)%index = eindex3D(ix+param%sx,iy+param%sy,iz+param%sz)
               gidvalue(ix)%value = real(permeab(0,ix,iy,iz),kind=4)
            end do
            call MPI_File_write( ounit, gidvalue, int(sizeof(gidvalue)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
         end do
      end do

      offset = ngrid_after * sizeof(gidvalue(1))
      call MPI_File_seek( ounit, offset, MPI_SEEK_CUR, mpi_ierr )

      call MPI_File_write_all( ounit, -1, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )   
      call write_string('End Values')

      call MPI_File_sync( ounit, mpi_ierr )

      return

      !--- Close file.
      entry doneResult
      
      call MPI_File_close( ounit, mpi_ierr )

   contains
      !>
      !! @brief output a string in GiD format.
      !! @param [in] str a string to be output.
      subroutine write_string( str )
         character(len=*), intent(in) :: str

         call MPI_File_write_all( ounit, len_trim(str)+1,    1,            MPI_INTEGER, MPI_STATUS_IGNORE, mpi_ierr )
         call MPI_File_write_all( ounit, trim(str)//char(0), len_trim(str)+1, MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
      end subroutine write_string

      !>
      !! @brief return serial index of a given element in global grid.
      !! @param [in] ix_g global index along x-axis
      !! @param [in] iy_g global index along y-axis
      !! @param [in] iz_g global index along z-axis
      !! @retval serial index.
      function eindex3D( ix_g, iy_g, iz_g ) result(index)
         integer, intent(in) :: ix_g, iy_g, iz_g
         integer :: index
         
         index = (ix_g-1) + (iy_g-1)*param%nx_g + (iz_g-1)*param%nx_g*param%ny_g + 1
      end function eindex3D

      !>
      !! @brief return serial index of a given mesh point in global grid.
      !! @param [in] ix_g global index along x-axis
      !! @param [in] iy_g global index along y-axis
      !! @param [in] iz_g global index along z-axis
      !! @retval serial index.
      function mindex3D( ix_g, iy_g, iz_g ) result(index)
         integer, intent(in) :: ix_g, iy_g, iz_g
         integer :: index

         index = (ix_g-1) + (iy_g-1)*(param%nx_g+1) + (iz_g-1)*(param%nx_g+1)*(param%ny_g+1) + 1 ! with extra point.
      end function mindex3D
   end subroutine saveResult

   !>
   !! @brief save fields in field file.
   !! @param [in] time current time.
   subroutine saveField( time )
      !--- Dummy arguments.
      real(8), intent(in) :: time
      !--- Local variables.
      integer :: ounit, status
      integer(8) :: offset

      type Header_type
         character(len=8) :: magic
         character(len=8) :: version
         real(8)    :: time
         integer(2) :: nx_g, ny_g, nz_g, dumm
      end type Header_type
      type(Header_type) :: head

      call MPI_File_delete( param%filename_last, MPI_INFO_NULL, mpi_ierr )
      call MPI_File_open( MPI_COMM_WORLD, param%filename_last, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, ounit, status )
      if( status/=0 ) then
         write(*,*) "Error: can not create field file: ", trim(param%filename_last)
         call exit(1)
      end if

      head%magic = "wPFAWHFG"
      head%version = "20250304"
      head%time  = time
      head%nx_g  = param%nx_g
      head%ny_g  = param%ny_g
      head%nz_g  = param%nz_g
      head%dumm  = 0

      offset = 0
      call MPI_File_write_all( ounit, head, int(sizeof(head)), MPI_BYTE, MPI_STATUS_IGNORE, mpi_ierr )
      offset = offset + sizeof(head)

      call MPI__File_write_body( ounit, offset, permeab_rock )
      offset = offset + param%ngrid_g*sizeof(permeab_rock(1,1,1))

      call MPI__File_write_body( ounit, offset, phase )
      offset = offset + param%ngrid_g*sizeof(phase(1,1,1))

      call MPI__File_write_body( ounit, offset, concent )
      offset = offset + param%ngrid_g*sizeof(concent(1,1,1))

      call MPI_File_close( ounit, mpi_ierr )
   end subroutine saveField

   !>
   !! @brief free arrays of fields.
   subroutine freeField
      deallocate( permeab, permeab_rock, pressure, velocity )
      deallocate( phase, concent, dphasedt, dconcentdt )
      deallocate( lw, well )
   end subroutine freeField

   function getPressure( ix_g, iy_g, iz_g ) result(p)
      integer, intent(in) :: ix_g, iy_g, iz_g
      real(8) :: p
      real(8) :: p_mpi

      if( param%sx+1<=ix_g .and. ix_g<=param%sx+param%nx .and. &
         param%sy+1<=iy_g .and. iy_g<=param%sy+param%ny .and. &
         param%sz+1<=iz_g .and. iz_g<=param%sz+param%nz ) then
         p = pressure(ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)
      else
         p = 0.0d0
      end if

      p_mpi = p
      call MPI_Allreduce( p_mpi, p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_ierr )
   end function getPressure

end module wPFAW3D

!>
!! @brief main routine of this program.
program main
   use wPFAW3D
   use, intrinsic :: ieee_arithmetic

   implicit none

   !---  Variables --------------------------------------------------------
   integer :: itime, i
   integer :: ix_g, iy_g, iz_g
   real(8) :: time, time_res, time_dat
   real(8) :: taumin
   character(len=256) :: filename, filename_dat, filename_res

   call MPI__Initialize
   call PetscInitialize( mpi_ierr )

   !--- Get filename of input or config file.
   if( command_argument_count()/=1 ) then
      write(*,*) "Error: no input file is given."
      call exit(1)
   end if
   call get_command_argument( 1, filename )

   !--- Load input or config file.
   i=index(filename,".inp")
   if( i>0 ) then
      call loadInput( filename )
      param%filename_last = filename(1:i-1)//".last.fld"
   else
      i=index(filename,".cfg")
      if( i>0 ) then
         call loadConfig( filename )
      else
         write(*,*) "Error: input file is neither inp or cfg."
         call exit(1)
      end if
   end if

   !--- Set filenames of data file and result file.
   filename_dat = filename(1:i-1)//".dat"
   filename_res = filename(1:i-1)//".res"

   !--- Start the time.
   itime    = 0                    ! current time step.
   time     = param%time_start     ! current time.
   time_res = time+param%dtime_res ! time to output at the next time.
   time_dat = time+param%dtime_dat ! time to output at the next time.

   !--- Write other information.
   call saveResult( filename_res, time )
   call showProgress( filename_dat, time, itime, getPressure(1,1,1) )

   !--- Main Calculation.
   do while( time<param%time_end )
      itime = itime+1 ! update step.
      if( itime>=3 ) then ! after 3rd step.
         param%dtime = param%tfac*0.5d0*param%dx**2/param%epsilon**2 &
            * min( param%epsilon**2/param%diffusion, taumin ) ! modify time step.
      end if
      time  = time+param%dtime ! update time.


      if (itime == 6) then
         open(5, file='phasecheck_ihou.out')
         write(5,*)"itime=",itime
         do iz_g = 1, 10
            write(5, *) iz_g
            do iy_g = 1, 10
               ! 1行に ix_g = 5〜14 の値を並べる
               write(5, '(11e12.4)') ( real(phase(ix_g, iy_g, iz_g)), ix_g = 1, 10 )
            end do
            write(5,*)
         enddo
      endif


      !--- Calculate time derivatives of fields.
      call calcPermeability
      call solvePressure
      call calcVelocity
      call calcPhaseDeriv( taumin )
      call calcConcentDeriv




      !    write(filename,"('field.midout.',i2.2)") mpi_rank
      !  open( unit=42, file=filename )   

      ! do ix_g=1, param%nx_g
      !    iy_g=ix_g
      !    iz_g=ix_g
      !    if( param%sx+1<=ix_g .and. ix_g<=param%sx+param%nx .and. &
      !       param%sy+1<=iy_g .and. iy_g<=param%sy+param%ny .and. &
      !       param%sz+1<=iz_g .and. iz_g<=param%sz+param%nz ) then
      !       write(42,"(i4,4e16.8)") ix_g, &
      !          real(pressure( ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
      !          real(permeab(0,ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
      !          real(phase(    ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
      !          real(concent(  ix_g-param%sx,iy_g-param%sy,iz_g-param%sz))
      !    end if
      ! end do
      ! write(42, *)



      !--- Evolve fields.
      call evolvePhase
      call evolveConcent

      !--- Print calculation progress.
      if( time>time_dat ) then ! time to output.
         call showProgress( filename_dat, time, itime, getPressure(param%ix_ref,param%iy_ref,param%iz_ref) )
         time_dat = time_dat+param%dtime_dat ! update time to output at the next time.
      end if

      !--- Print calculation results.
      if( time>time_res ) then ! time to output.
         call saveResult( filename_res, time )
         time_res = time_res+param%dtime_res ! update time to output at the next time.
      end if
      ! write(*,*)"itime = ",itime
   end do

   !--- For debug
   write(filename,"('field.out.',i2.2)") mpi_rank
   open( unit=100, file=filename )

   do ix_g=1, param%nx_g
      iy_g=ix_g
      iz_g=ix_g
      if( param%sx+1<=ix_g .and. ix_g<=param%sx+param%nx .and. &
         param%sy+1<=iy_g .and. iy_g<=param%sy+param%ny .and. &
         param%sz+1<=iz_g .and. iz_g<=param%sz+param%nz ) then
         write(100,"(i4,4e16.8)") ix_g, &
            real(pressure( ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
            real(permeab(0,ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
            real(phase(    ix_g-param%sx,iy_g-param%sy,iz_g-param%sz)), &
            real(concent(  ix_g-param%sx,iy_g-param%sy,iz_g-param%sz))
      end if
   end do

   close(100)
   !--- For debug

   !--- Finish this program.
   call doneProgress
   call doneResult
   call saveField( time )
   call freeField

   call PetscFinalize( mpi_ierr )
   call MPI__Finalize
end program main
