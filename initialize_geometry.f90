subroutine initialize_geometry (geometry,params)

  ! Subroutine to initialize the gemetry of the problem
  ! such as node locations, velocities, erosional properties, etc

  use definitions

  implicit none

  type (parm) params
  type (geom) geometry
  integer nx,ny,n
  integer, dimension(:), allocatable :: seed
  double precision xl,yl,zl,xx, spacing
  integer ij,i,j,fix


  call time_in ('initialize_geometry')

  ! In this beta version we will position the nodes on a quasi regular rectangular grid
  ! of size xl by yl and nx by ny nodes; the location of the nodses will be randomly "shaken"
  ! zl is the maximum height of the initial landscape


nx = geometry%nx
ny = geometry%ny
xl =  geometry%xl
yl =  geometry%yl
zl =  geometry%zl

  geometry%nnode=nx*ny ! computes total number of nodes
  geometry%nnode_max=geometry%nnode*5
  !geometry%nnmax=20 ! nnmax is the maximum number of neighbour nodes per node; 12 is usually enough
  geometry%ndivide_max = geometry%nnode_max*(geometry%nnmax-1) !max number of divides
  geometry%ndivide=0 ! don't change this
  geometry%ncapture=0 ! don't cahnge this


  print*, ' '
  print*, 'Grid parameters: '
  print*, 'nx = ',nx
  print*, 'ny = ',ny
  print*, 'xl = ',xl
  print*, 'yl = ',yl
  print*, 'zl = ',zl
  print*, ' '

  allocate (geometry%x(geometry%nnode_max),geometry%y(geometry%nnode_max),geometry%z(geometry%nnode_max))
  allocate (geometry%xdiv(geometry%ndivide_max),geometry%ydiv(geometry%ndivide_max),&
       &geometry%zdiv(geometry%ndivide_max))
  allocate (geometry%u(geometry%nnode_max),geometry%v(geometry%nnode_max),geometry%w(geometry%nnode_max))
  allocate (geometry%fix(geometry%nnode_max),geometry%surface(geometry%nnode_max))
  allocate (geometry%discharge(geometry%nnode_max),geometry%precipitation(geometry%nnode_max))
  allocate (geometry%k(geometry%nnode_max),geometry%catchment(geometry%nnode_max))
  allocate (geometry%nb(geometry%nnode_max),geometry%erosion_rate(geometry%nnode_max))
  allocate (geometry%nt(geometry%nnode_max))
  allocate (geometry%nn(geometry%nnmax,geometry%nnode_max))
  allocate (geometry%nnodetri(geometry%nnmax,geometry%nnode_max))
  allocate (geometry%nndivtri(geometry%ndivide_max,2))
  allocate (geometry%nndivnode(geometry%ndivide_max,2))
  allocate (geometry%sediment_flux(geometry%nnode_max), geometry%chi(geometry%nnode_max))
  allocate (geometry%strahler(geometry%nnode_max))
  allocate (geometry%surface_share(geometry%nnmax,geometry%nnode_max))
  if (params%transient_divide) then
     allocate(geometry%erosion_rate_history(geometry%nnode_max,params%num_bins))
  endif

  geometry%nb=0 ! don't change this
  geometry%nt=0

  geometry%k=params%k_scalar1    ! fluvial erosion constant
  geometry%erosion_rate=0.1d-4 ! sets erosion rate to 0
  geometry%sediment_flux=0.d0 ! sets sediment flux to 0
  geometry%discharge=0.d0
  geometry%precipitation=0.d0
  geometry%catchment=0.d0
  geometry%fix=0.d0
  geometry%surface=0.d0
  geometry%strahler=0.d0
  geometry%surface_share=0.d0


  ! generate random numbers in x, y and z
  call random_seed(SIZE = n)
  allocate (seed(n))
  seed=6
  call random_seed(PUT= seed(1:n))
  call random_number (geometry%z)
  call random_number (geometry%x)
  call random_number (geometry%y)
  geometry%x=geometry%x*2.d0-1.d0
  geometry%y=geometry%y*2.d0-1.d0
  geometry%z=.01+geometry%z*zl

  ! generates the grid

  !call random_number(xx)
  ij=0
  do j=1,ny
     do i=1,nx
        ij=ij+1
        fix=0
        if (i.eq.1 .or. i.eq.nx .or. j.eq.1 .or. j.eq.ny) fix=1
        geometry%x(ij)=(geometry%x(ij)/(nx-1)/4.*(1-fix)+float(i-1)/(nx-1))*xl
        geometry%y(ij)=(geometry%y(ij)/(ny-1)/4.*(1-fix)+float(j-1)/(ny-1))*yl
        geometry%fix(ij)=0 ! geometry%fix=1 means that the node is at a fixed base level
        !   call random_number(xx)
        !    xx=0.1
        !call random_number (xx)
        !if (i.eq.1 .or. i.eq.nx .or. j.eq.1 .or. j.eq.ny) geometry%fix(ij)=1
        if(i.eq.1)geometry%fix(ij) = 1
        if(i.eq.nx)geometry%fix(ij) = 2
        if(j.eq.1)geometry%fix(ij) = 3
        if(j.eq.ny)geometry%fix(ij) = 4


        geometry%u(ij)=0.d0
        geometry%v(ij)=0.d0
        if (geometry%fix(ij).ge.1) then
           geometry%z(ij)=0.d0
           geometry%erosion_rate(ij)=0.d0*params%uplift_scalar2
        endif
     enddo
  enddo




  if (params%transient_divide) then
     geometry%erosion_rate_history=0.d0
  endif



   ! housekeeping - initialize all arrayse
   geometry%precipitation=params%rainfall_height
    geometry%u=0
    geometry%v=0
    geometry%w=params%uplift_scalar1
   geometry%chi=0






  call time_out ('initialize_geometry')

  return

end subroutine initialize_geometry
