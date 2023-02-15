subroutine find_hack (geometry,network,stack)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  integer i,j,k
  double precision l
  double precision,dimension(:),allocatable::contributing_area, max_length, horton_length

  call time_in ('find_hack')

  allocate(contributing_area(geometry%nnode), max_length(geometry%nnode), horton_length(geometry%nnode))

  contributing_area=geometry%surface
  max_length=0.
  do i=stack%nnode,1,-1
     j=stack%order(i)
     k=network%receiver(j)
     if (k.ne.0) then
        contributing_area(k) = contributing_area(k) + contributing_area(j)
        l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
        max_length(k) = max(max_length(k),max_length(j)+l)
     endif
  enddo
  
  open (40,file='HackData',status='unknown')
  do i=1,geometry%nnode
     write(40,'(2i7,2f16.3,i3)') i, geometry%catchment(i), contributing_area(i), max_length(i), geometry%strahler(i)
  enddo
  close (40,err=1232)
  open (41,file='HortonData',status='unknown')
  horton_length = 0.
  do i=stack%nnode,1,-1
      j=stack%order(i)
      if (network%receiver(j) .eq. 0) then
         if (horton_length(j).gt.0.d0) then
            write(41,'(2i7,2f16.3,i3)') j, geometry%catchment(j), contributing_area(j), horton_length(j), geometry%strahler(j)
         endif
      else
         k=network%receiver(j)
         l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
         if (geometry%strahler(j).eq.geometry%strahler(k)) then    
            horton_length(k) = horton_length(j) + l
         else
            horton_length(j) =  horton_length(j) + l
            write(41,'(2i7,2f16.3,i3)') j, geometry%catchment(j), contributing_area(j), horton_length(j), geometry%strahler(j)
         endif
      endif
   enddo
  close (41,err=1232)
  deallocate(contributing_area, max_length,horton_length)
  call time_out ('find_hack')

  return
1232 STOP 'problem closing HackData'
end subroutine find_hack
