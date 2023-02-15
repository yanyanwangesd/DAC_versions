subroutine write_ascii (geometry,params,network)
  use definitions

  implicit none

  type (geom) geometry
  type (netw) network 
  type (parm) params
  integer i,icount,ii
  character cs*8
  


  icount=params%istep
 

  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'

  open (70,file='ASCII/coordinates'//cs,status='unknown')
  open (71,file='ASCII/erosion_rate'//cs,status='unknown')
  open (72,file='ASCII/sediment_flux'//cs,status='unknown')
  open (73,file='ASCII/catchment'//cs,status='unknown')
  open (74,file='ASCII/river_network'//cs,status='unknown')

  
  do i=1,geometry%nnode
     write(70,'(3f16.3)') geometry%x(i),geometry%y(i),geometry%z(i)
  enddo
  
  do i=1,geometry%nnode
     write(71,'(e12.4)') geometry%erosion_rate(i)
  enddo

  do i=1,geometry%nnode
     write(72,'(e12.4)') geometry%sediment_flux(i)
  enddo

  do i=1,geometry%nnode
     write(73,'(i5)') geometry%catchment(i)
  enddo
  
  do i=1,geometry%nnode
     ii=network%receiver(i)
     if (network%receiver(i).eq.0) ii=i
     write(74,'(9I10)') i,ii
  enddo

  close (70,err=1231)
  close (71,err=1231)
  close (72,err=1231)
  close (73,err=1231)
  close (74,err=1231)
  return
1231 STOP 'problem closing file in VTK'
end
