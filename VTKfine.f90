subroutine VTKfine (geometry,delaunay,params,network,stack,landscape)
  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay
  type (netw) network
  type (stck) stack
  type (parm) params
  TYPE (lands) landscape
  integer i,icount,iordermin,nlinks,ii,iorder, j, k
  character cs*8
  double precision vex,v1,v2
  double precision max_erosion_rate,  min_precipitation


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



  open (35,file='RUN5/LandscapeFine'//cs//'.vtk',status='unknown')


  write(35,'(a)')'# vtk DataFile Version 3.0'
  write(35,'(a)')'Landscape'
  write(35,'(a)')'ASCII'
  write(35,'(a)')'DATASET UNSTRUCTURED_GRID'

  write(35,'(a7,i10,a6)')'POINTS ',landscape%nfinenode,' float'
  do i=1,landscape%nfinenode
     write(35,'(3f16.3)') landscape%xx(i),landscape%yy(i),landscape%zz(i)
  enddo

  write(35,'(A6, 2I10)') 'CELLS ',landscape%nfinetri,4*landscape%nfinetri
  do i=1,landscape%nfinetri
     write(35,'(9I10)')3,landscape%fineicon(1:3,i)-1
  enddo

  write(35,'(A11, I10)') 'CELL_TYPES ',landscape%nfinetri
  do i=1,landscape%nfinetri
     write(35,'(I2)')5 ! octree  (8 nodes)
  enddo

  write(35,'(a11,i10)')'POINT_DATA ',landscape%nfinenode
  write(35,'(a)')'SCALARS height float 1'
  write (35,'(a)')'LOOKUP_TABLE default'
  do i=1,landscape%nfinenode
     write(35,'(F8.1)') landscape%zz(i)
  enddo

  write(35,'(a)')'SCALARS erosion_rate float 1'
  write (35,'(a)')'LOOKUP_TABLE default'
  do i=1,landscape%nfinenode
     write(35,'(e12.4)') landscape%edot(i)
  enddo


  close (35,err=1339)






  return
1339 STOP 'problem closing file in VTKfine'
end subroutine VTKfine

