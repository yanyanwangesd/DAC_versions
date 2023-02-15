subroutine reinit_geometry (geometry,params,network,delaunay)

  ! Subroutine to initialize the gemetry of the problem
  ! such as node locations, velocities, erosional properties, etc

  use definitions

  implicit none

  type (parm) params
  type (geom) geometry
  type (netw) network
  type (del) delaunay
  integer i,j, nnode,xmin,xmax,ymin,ymax
  character*6   cs

  call time_in ('reinit_geometry')

  !read (5,'(a)') params%title
  !print*, 'reading restart file ',params%title, ' length is ', len(params%title)
  !print*, 'file number is ', title(8:15), 'length is ',len(title(8:15)) 
  !print*, 'trimmed file number is ', trim(title(8:15)), 'trimmed length is ',len(trim(title(8:15))) 
  !read(title(8:11),'(i4)') params%num_restart
  !params%num_restart = params%num_restart + 1
  !print*, params%num_restart
  !print*, params%num_restart + 1
  !print*, '******'
  !open (5,file='GeoFile',status='unknown')

  write(cs,'(I6)') params%num_restart
  if (params%num_restart.lt.10)        cs(1:5)='00000'
  if (params%num_restart.lt.100)       cs(1:4)='0000'
  if (params%num_restart.lt.1000)      cs(1:3)='000'
  if (params%num_restart.lt.10000)     cs(1:2)='00'
  if (params%num_restart.lt.100000)    cs(1:1)='0'

  
  open(unit=51,file= 'RESTART/GeoFile'//cs,form='unformatted')
  print*, 'reading restart file ', 'RESTART/GeoFile'//cs
  !params%num_restart = params%num_restart + 1

  read(51) params%time
  read(51) params%istep
  read(51) nnode
  

  params%tfinal = params%time+params%tfinal
  print*, 'endtime = ',params%tfinal
  geometry%nnode=nnode ! computes total number of nodes
  geometry%nnode_max=geometry%nnode*51
  !geometry%nnmax=12 ! nnmax is the maximum number of neighbour nodes per node; 12 is usually enough
  geometry%ndivide_max = geometry%nnode_max*(geometry%nnmax-1) !max number of divides
  geometry%ndivide=0 ! don't change this
  geometry%ncapture=0 ! don't cahnge this

  allocate (geometry%x(geometry%nnode_max),geometry%y(geometry%nnode_max),geometry%z(geometry%nnode_max))
  do i =1,nnode
     read(51) geometry%x(i),geometry%y(i),geometry%z(i)
!     print*, geometry%x(i),geometry%y(i),geometry%z(i)
  enddo




  xmin = minval(geometry%x)
  xmax = maxval(geometry%x)
  ymin = minval(geometry%y)
  ymax = maxval(geometry%y)





  allocate (geometry%xdiv(geometry%ndivide_max),geometry%ydiv(geometry%ndivide_max),&
       &geometry%zdiv(geometry%ndivide_max))	
  allocate (geometry%u(geometry%nnode_max),geometry%v(geometry%nnode_max),geometry%w(geometry%nnode_max))
  allocate (geometry%fix(geometry%nnode_max),geometry%surface(geometry%nnode_max))
  allocate (geometry%discharge(geometry%nnode_max),geometry%precipitation(geometry%nnode_max))
  allocate (geometry%k(geometry%nnode_max),geometry%catchment(geometry%nnode_max))
  allocate (geometry%chi(geometry%nnode_max))
  allocate (geometry%nb(geometry%nnode_max),geometry%erosion_rate(geometry%nnode_max))
  allocate (geometry%nt(geometry%nnode_max))
  allocate (geometry%nn(geometry%nnmax,geometry%nnode_max))
  allocate (geometry%nnodetri(geometry%nnmax,geometry%nnode_max))
  allocate (geometry%nndivtri(geometry%ndivide_max,2))
  allocate (geometry%nndivnode(geometry%ndivide_max,2))
  allocate (geometry%sediment_flux(geometry%nnode_max))
  allocate (geometry%strahler(geometry%nnode_max))
  allocate (geometry%surface_share(geometry%nnmax,geometry%nnode_max))
  if (params%transient_divide) then
     allocate(geometry%erosion_rate_history(geometry%nnode_max,params%num_bins))
  endif
  
  if (params%transient_divide) then
    do i = 1,geometry%nnode
       do j = 1,params%num_bins
          read(51) geometry%erosion_rate_history(i,j)
       enddo
    enddo
 endif
  
 geometry%nb=0 ! don't change this
 geometry%nt=0



 do i=1,geometry%nnode
    read(51) geometry%erosion_rate(i),geometry%sediment_flux(i)
 enddo
  

  geometry%k=params%k_scalar1    ! fluvial erosion constant
 


  do i = 1,nnode
     if (geometry%x(i).eq.xmin.or.geometry%x(i).eq.xmax.or.geometry%y(i).eq.ymin.or.geometry%y(i).eq.ymax) then
        geometry%fix(i) = 1
     else
        geometry%fix(i) = 0
     endif
  enddo


  

  read(51) geometry%ndivide
  do i = 1,geometry%ndivide
     read(51) geometry%xdiv(i), geometry%ydiv(i), geometry%zdiv(i)
  enddo
  do i = 1,geometry%ndivide
     read(51) geometry%nndivnode(i,1), geometry%nndivnode(i,2), geometry%nndivtri(i,1),geometry%nndivtri(i,2)
  enddo


  network%nnode=geometry%nnode
  network%ndonmax=geometry%nnmax
  allocate (network%donors(network%ndonmax,geometry%nnode_max))
  allocate (network%receiver(geometry%nnode_max))
  allocate (network%lakes_catch(geometry%nnode_max))
  allocate (network%lakes(geometry%nnode_max))
  allocate (network%ndon(geometry%nnode_max))

  do i=1,network%nnode
     read(51) network%ndon(i)
     if (network%ndon(i).ne.0) then
        read(51) (network%donors(j,i),j=1,network%ndon(i))
     endif
  enddo
  do i=1,network%nnode
     read(51) network%receiver(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%discharge(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%erosion_rate(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%sediment_flux(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%catchment(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%chi(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%precipitation(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%surface(i)
  enddo
  do i=1,geometry%nnode
     read(51) geometry%fix(i)
  enddo
 


  read(51) delaunay%ntriangles

  

  allocate (delaunay%icon(3,delaunay%ntriangles),delaunay%neighbours(3,delaunay%ntriangles),delaunay%numdivides(4,delaunay%ntriangles),&
       &delaunay%centers(3,delaunay%ntriangles))

  do i=1,delaunay%ntriangles
     read(51) delaunay%icon(1,i), delaunay%icon(2,i), delaunay%icon(3,i)
  enddo

  do i=1,delaunay%ntriangles
     read(51) delaunay%neighbours(1,i), delaunay%neighbours(2,i), delaunay%neighbours(3,i)
  enddo
  do i=1,delaunay%ntriangles
     read(51) delaunay%numdivides(1,i), delaunay%numdivides(2,i), delaunay%numdivides(3,i),delaunay%numdivides(4,i)
  enddo
  do i=1,delaunay%ntriangles
     read(51) delaunay%centers(1,i), delaunay%centers(2,i), delaunay%centers(3,i)
  enddo



   ! housekeeping - initialize all arrayse
   geometry%precipitation=params%rainfall_height
    geometry%u=0
    geometry%v=0
    geometry%w=params%uplift_scalar1



  !do i = 1,nnode
  !   read(51) geometry%nb(i)
  !   read(51) (geometry%nn(j,i),j=1,geometry%nb(i))
  !enddo

  close (51,err=1233)



  call time_out ('reinit_geometry')

  return

1233 STOP 'problem closing GeoFile'
end subroutine reinit_geometry
