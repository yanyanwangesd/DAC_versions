subroutine VTK (geometry,delaunay,params,network,stack)
  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay
  !type (del) rddelaunay
  type (netw) network
  type (stck) stack
  type (parm) params
  integer i,icount,iordermin,nlinks,ii,iorder, j, k
  character cs*8
  double precision vex,v1,v2
  double precision max_erosion_rate,  min_precipitation

!!$  max_erosion_rate = 0.0
!!$  min_precipitation = params%rainfall_height
!!$  do i = 1,geometry%nnode
!!$     max_erosion_rate = max(max_erosion_rate,geometry%erosion_rate(i))
!!$     min_precipitation = min( min_precipitation, geometry%precipitation(i))
!!$  enddo
!!$  print*, 'max erosion rate is', max_erosion_rate
!!$  print*, 'minprecipitation is',  min_precipitation 
  
  v1=minval(geometry%z)
  v2=maxval(geometry%z)
  !  print*,'Topo min-max in VTK',v1,v2
  print*,'time=',params%time/1.e3,'kyrs'
  !  print*,'cumulative captures',geometry%ncapture

  icount=params%istep
  vex=1.d0

  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'

  open (30,file='RUN5/Landscape'//cs//'.vtk',status='unknown')
  open (31,file='RUN5/rivers'//cs//'.vtk',status='unknown')
  open (32,file='RUN5/divides'//cs//'.vtk',status='unknown')
  !open (33,file='RUN5/RDLandscape'//cs//'.vtk',status='unknown')

  write(30,'(a)')'# vtk DataFile Version 3.0'
  write(30,'(a)')'Landscape'
  write(30,'(a)')'ASCII'
  write(30,'(a)')'DATASET UNSTRUCTURED_GRID'

  write(30,'(a7,i10,a6)')'POINTS ',geometry%nnode,' float'
  do i=1,geometry%nnode
     write(30,'(3f16.3)') geometry%x(i),geometry%y(i),geometry%z(i)*vex
  enddo

  write(30,'(A6, 2I10)') 'CELLS ',delaunay%ntriangles,4*delaunay%ntriangles
  do i=1,delaunay%ntriangles
     write(30,'(9I10)')3,delaunay%icon(1:3,i)-1
  enddo

  write(30,'(A11, I10)') 'CELL_TYPES ',delaunay%ntriangles
  do i=1,delaunay%ntriangles
     write(30,'(I2)')5 ! octree  (8 nodes)
  enddo

  write(30,'(a11,i10)')'POINT_DATA ',geometry%nnode
  write(30,'(a)')'SCALARS height float 1'
  write (30,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%nnode
     write(30,'(F8.1)') geometry%z(i)
  enddo

  write(30,'(a)')'SCALARS discharge float 1'
  write (30,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%nnode
     write(30,'(e12.4)') geometry%discharge(i)
  enddo

  write(30,'(a)')'SCALARS erosion_rate float 1'
  write (30,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%nnode
     write(30,'(e12.4)') geometry%erosion_rate(i)
  enddo

  write(30,'(a)')'SCALARS sediment_flux float 1'
  write (30,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%nnode
     write(30,'(e12.4)') geometry%sediment_flux(i)
  enddo

  write(30,'(a)')'SCALARS catchment float 1'
  write (30,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%nnode
     write(30,'(i5)') geometry%catchment(i)
  enddo
  !pause
  close (30,err=1231)

  ! Writing vtk files with the rivers
  iordermin=1
  nlinks=geometry%nnode

  write(31,'(a)')'# vtk DataFile Version 3.0'
  write(31,'(a)')'rivers'
  write(31,'(a)')'ASCII'
  write(31,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(31,'(a7,i10,a6)')'POINTS ',geometry%nnode,' float'
  do i=1,geometry%nnode
     write(31,'(3f16.3)') geometry%x(i),geometry%y(i),geometry%z(i)!+geometry%z(i)/10.
  enddo

  write(31,'(a6, 2I10)') 'CELLS ',nlinks,(1+2)*nlinks
  do i=1,geometry%nnode
     ii=network%receiver(i)
     if (network%receiver(i).eq.0) ii=i
     write(31,'(9I10)') 2,i-1,ii-1
  enddo

  write(31,'(A11, I10)') 'CELL_TYPES ',nlinks
  do i=1,nlinks
     write(31,'(I2)') 3 ! octree  (2 nodes)
  end do

  write(31,'(a11,i10)')'POINT_DATA ',geometry%nnode
  write(31,'(a)')'SCALARS discharge float 1'
  write (31,'(a)')'LOOKUP_TABLE default'

  do i=1,geometry%nnode
     write(31,'(e12.4)') geometry%discharge(i)
  enddo


  write(31,'(a)')'SCALARS chi float 1'
  write (31,'(a)')'LOOKUP_TABLE default'

  do i=1,geometry%nnode
     write(31,'(e12.4)') geometry%chi(i)
  enddo


!!$  write(31,'(a)')'SCALARS strahler float 1'
!!$  write (31,'(a)')'LOOKUP_TABLE default'
!!$
!!$  do i=1,geometry%nnode
!!$     write(31,'(i3)') geometry%strahler(i)
!!$  enddo


  write(31,'(a11,i10)')'CELL_DATA ',nlinks
  write(31,'(a)')'SCALARS strahler float 1'
  write (31,'(a)')'LOOKUP_TABLE default'

  do i=1,geometry%nnode
     !ii=network%receiver(i)
     !if (network%receiver(i).eq.0) then 
     write(31,'(i3)') geometry%strahler(i)
     !else 
     !   write(31,'(i3)') geometry%strahler(ii)
     !endif
  enddo

 
  close (31,err=1231)


  ! Writing vtk files with the divides
  print*, 'num of divides is ',geometry%ndivide, ' num of triangles is ',delaunay%ntriangles
  nlinks = 0
  do i = 1,delaunay%ntriangles
     if (delaunay%numdivides(1,i).eq.1) then
        nlinks=nlinks + 1
     elseif (delaunay%numdivides(1,i).eq.2) then 
        !nlinks=nlinks + 1
        nlinks=nlinks + 2
     elseif (delaunay%numdivides(1,i).eq.3) then
        nlinks=nlinks + 3
     else
!        Print*, 'triangle',i,'is having',delaunay%numdivides(1,i),'divides'
     endif
  enddo

  write(32,'(a)')'# vtk DataFile Version 3.0'
  write(32,'(a)')'divides'
  write(32,'(a)')'ASCII'
  write(32,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(32,'(a7,i10,a6)')'POINTS ',geometry%ndivide+delaunay%ntriangles,' float'
  do i=1,geometry%ndivide
     write(32,'(3f16.3)') geometry%xdiv(i),geometry%ydiv(i),geometry%zdiv(i)*vex
  enddo
  do i=1,delaunay%ntriangles
     write(32,'(3f16.3)') delaunay%centers(1,i),delaunay%centers(2,i),delaunay%centers(3,i)
  enddo


  write(32,'(a6, 2I10)') 'CELLS ',nlinks,(1+2)*(nlinks)
  do i=1,delaunay%ntriangles
     j=delaunay%numdivides(1,i)
     if (j.eq.1) then
        write(32,'(9I10)') 2, geometry%ndivide + i - 1, delaunay%numdivides(2,i)-1 ! connect centroid to the divide
     elseif (j.eq.2) then
        !write(32,'(9I10)') 2, delaunay%numdivides(2,i)-1,delaunay%numdivides(3,i)-1 ! connect the two divides
        do k=1,2
           write(32,'(9I10)') 2, geometry%ndivide + i - 1, delaunay%numdivides(k+1,i)-1 !connect centroid to the three divides
        enddo
     elseif (j.eq.3) then
        do k=1,3
           write(32,'(9I10)') 2, geometry%ndivide + i - 1, delaunay%numdivides(k+1,i)-1 !connect centroid to the three divides
        enddo
     endif
  enddo

 

  write(32,'(A11, I10)') 'CELL_TYPES ',nlinks
  do i=1,nlinks
     write(32,'(I2)') 3 ! two nodes
  end do

  write(32,'(a11,i10)')'POINT_DATA ',geometry%ndivide+delaunay%ntriangles
  write(32,'(a)')'SCALARS height float 1'
  write (32,'(a)')'LOOKUP_TABLE default'
  do i=1,geometry%ndivide
     write(32,'(F8.1)') geometry%zdiv(i)
  enddo
  do i=1,delaunay%ntriangles
     write(32,'(F8.1)') delaunay%centers(3,i) 
  enddo

  close (32,err=1231)



  !writing vtk files for the landscape of rivers and divides
!!$
!!$
!!$  write(33,'(a)')'# vtk DataFile Version 3.0'
!!$  write(33,'(a)')'RDLandscape'
!!$  write(33,'(a)')'ASCII'
!!$  write(33,'(a)')'DATASET UNSTRUCTURED_GRID'
!!$
!!$  if (rddelaunay%ntriangles.ne.0) then
!!$
!!$     write(33,'(a7,i10,a6)')'POINTS ',geometry%nnode+geometry%ndivide,' float'
!!$     do i=1,geometry%nnode
!!$        write(33,'(3f16.3)') geometry%x(i),geometry%y(i),geometry%z(i)*vex
!!$     enddo
!!$     do i=1,geometry%ndivide
!!$        write(33,'(3f16.3)') geometry%xdiv(i),geometry%ydiv(i),geometry%zdiv(i)*vex
!!$     enddo
!!$
!!$     write(33,'(A6, 2I10)') 'CELLS ',rddelaunay%ntriangles,4*rddelaunay%ntriangles
!!$     do i=1,rddelaunay%ntriangles
!!$        write(33,'(9I10)')3,rddelaunay%icon(1:3,i)-1
!!$     enddo
!!$
!!$     write(33,'(A11, I10)') 'CELL_TYPES ',rddelaunay%ntriangles
!!$     do i=1,rddelaunay%ntriangles
!!$        write(33,'(I2)')5 ! octree  (8 nodes)
!!$     enddo
!!$
!!$     write(33,'(a11,i10)')'POINT_DATA ',geometry%nnode+geometry%ndivide
!!$     write(33,'(a)')'SCALARS height float 1'
!!$     write (33,'(a)')'LOOKUP_TABLE default'
!!$     do i=1,geometry%nnode
!!$        write(33,'(F8.1)') geometry%z(i)
!!$     enddo
!!$     do i=1,geometry%ndivide
!!$        write(33,'(F8.1)') geometry%zdiv(i)
!!$     enddo
!!$  else
!!$     write(33,'(a7,i10,a6)')'POINTS ',geometry%nnode,' float'
!!$     do i=1,geometry%nnode
!!$        write(33,'(3f16.3)') geometry%x(i),geometry%y(i),geometry%z(i)*vex
!!$     enddo
!!$
!!$     write(33,'(A6, 2I10)') 'CELLS ',delaunay%ntriangles,4*delaunay%ntriangles
!!$     do i=1,delaunay%ntriangles
!!$        write(33,'(9I10)')3,delaunay%icon(1:3,i)-1
!!$     enddo
!!$
!!$     write(33,'(A11, I10)') 'CELL_TYPES ',delaunay%ntriangles
!!$     do i=1,delaunay%ntriangles
!!$        write(33,'(I2)')5 ! octree  (8 nodes)
!!$     enddo
!!$
!!$     write(33,'(a11,i10)')'POINT_DATA ',geometry%nnode
!!$     write(33,'(a)')'SCALARS height float 1'
!!$     write (33,'(a)')'LOOKUP_TABLE default'
!!$     do i=1,geometry%nnode
!!$        write(33,'(F8.1)') geometry%z(i)
!!$     enddo
!!$
!!$  endif
!!$
!!$  !pause
!!$  close (33,err=1231)
!!$


  return
1231 STOP 'problem closing file in VTK'
end
