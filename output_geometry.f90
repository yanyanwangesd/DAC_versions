subroutine output_geometry (geometry,params,network,delaunay)


  use definitions

  implicit none

  type (geom) geometry
  type (parm) params
  type (netw) network
  type (del) delaunay
  integer i, j
  character*4  cs


  call time_in ('output_geometry')
  
  write(cs,'(I4)') params%num_restart
  if (params%num_restart.lt.10)        cs(1:3)='000'
  if (params%num_restart.lt.100)       cs(1:2)='00'
  if (params%num_restart.lt.1000)      cs(1:1)='0'
  !if (params%num_restart.lt.10000)     cs(1:4)=''
  
  
  params%num_restart = params%num_restart + 1
  print*, 'writing binary file ', 'RESTART/GeoFile'//cs
  !open (50,file='RESTART/GeoFile'//cs,status='unknown')
  !open (50,file='GeoFile',status='unknown')
  open (50,file='RESTART/GeoFile'//cs,form = 'unformatted')
  write(50) params%time
  write(50) params%istep
  write(50) geometry%nnode
  do i=1,geometry%nnode
     write(50) geometry%x(i),geometry%y(i),geometry%z(i)
  enddo
  
 
 if (params%transient_divide) then
    do i = 1,geometry%nnode
       do j = 1,params%num_bins
          write(50) geometry%erosion_rate_history(i,j)
       enddo
    enddo
 endif
   
 
 do i=1,geometry%nnode
    write(50) geometry%erosion_rate(i),geometry%sediment_flux(i)
 enddo



  write(50) geometry%ndivide
  do i = 1,geometry%ndivide
     write(50) geometry%xdiv(i), geometry%ydiv(i), geometry%zdiv(i)
  enddo
  do i = 1,geometry%ndivide
     write(50) geometry%nndivnode(i,1), geometry%nndivnode(i,2), geometry%nndivtri(i,1),geometry%nndivtri(i,2)
  enddo


  do i=1,network%nnode
     write(50) network%ndon(i)
     if (network%ndon(i).ne.0) then
        write(50) (network%donors(j,i),j=1,network%ndon(i))
     endif
  enddo
  do i=1,network%nnode
     write(50) network%receiver(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%discharge(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%erosion_rate(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%sediment_flux(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%catchment(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%chi(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%precipitation(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%surface(i)
  enddo
  do i=1,geometry%nnode
     write(50) geometry%fix(i)
  enddo
  


  write(50) delaunay%ntriangles
  
  do i=1,delaunay%ntriangles
     write(50) delaunay%icon(1,i), delaunay%icon(2,i), delaunay%icon(3,i)
  enddo
  do i=1,delaunay%ntriangles
     write(50) delaunay%neighbours(1,i), delaunay%neighbours(2,i), delaunay%neighbours(3,i)
  enddo
  do i=1,delaunay%ntriangles
     write(50) delaunay%numdivides(1,i), delaunay%numdivides(2,i), delaunay%numdivides(3,i),delaunay%numdivides(4,i)
  enddo
  do i=1,delaunay%ntriangles
     write(50) delaunay%centers(1,i), delaunay%centers(2,i), delaunay%centers(3,i)
  enddo

  

  !   do i = 1,geometry%nnode
  !      write(50,*) geometry%nb(i)
  !      write(50,*) (geometry%nn(j,i),j=1,geometry%nb(i))
  !   enddo


  close (50,err=1233)


  call time_out ('output_geometry')

  return
1233 STOP 'problem closing GeoFile'
end subroutine output_geometry
