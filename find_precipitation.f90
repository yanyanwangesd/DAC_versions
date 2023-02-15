subroutine find_precipitation (geometry,params)

! computes orographic precipitation (not fully implemented yet)
! in this version precipitation is set to 1 (m/yr)

use definitions

implicit none

type (geom) geometry
type (parm) params
type (netw) network
type (stck) stack
double precision slopemax,slope,rainfall
double precision deltax,deltay,nx,ny,direction,yl
integer i,j,k

call time_in ('find_precipitation')
yl=maxval(geometry%y)-minval(geometry%y)

geometry%precipitation=params%rainfall_height


! hard code change in precip - 10-2010
!if(params%time.gt.params%tfinal/2.) geometry%precipitation=params%rainfall_height/10.


!geometry%precipitation=params%rainfall_height*((geometry%y-minval(geometry%y))/yl)

call time_out ('find_precipitation')
return











direction=30.d0
direction=direction/45.d0*atan(1.d0)
nx=cos(direction)
ny=sin(direction)

network%nnode=geometry%nnode
network%ndonmax=12
allocate (network%donors(network%ndonmax,network%nnode))
allocate (network%receiver(network%nnode))
allocate (network%ndon(network%nnode))

network%receiver=0
  do i=1,geometry%nnode
  if (geometry%nb(i).eq.0) stop 'Disconnected node in create_network'
!    if (geometry%fix(i).eq.0) then
    slopemax=-2.d0
      do j=1,geometry%nb(i)
      k=geometry%nn(j,i)
      deltax=geometry%x(i)-geometry%x(k)
      deltay=geometry%y(i)-geometry%y(k)
      slope=sqrt(deltax**2+deltay**2)
      slope=(deltax*nx+deltay*ny)/slope
        if (slope.gt.slopemax) then
        slopemax=slope
        network%receiver(i)=k
        endif
      enddo
!    endif
  enddo

  do i=1,geometry%nnode
  if (geometry%fix(i).ge.1 .and. geometry%fix(network%receiver(i)).ge.1) network%receiver(i)=0
  enddo

print*,'total number of nodes', geometry%nnode
print*,'number of nodes with a receiver ',count(network%receiver.ne.0)

network%donors=0
network%ndon=0
  do i=1,network%nnode
  k=network%receiver(i)
    if (k.ne.0) then
    network%ndon(k)=network%ndon(k)+1
    network%donors(network%ndon(k),k)=i
    endif
  enddo

stack%nnode=0
!print*,'number of donors ',sum(network%ndon)
!print*,count(network%ndon==0)
!call find_order (geometry,network,stack,delaunay,params)

geometry%precipitation=0.d0
  do i=1,stack%nnode
  j=stack%order(i)
  if (network%receiver(j).eq.0) rainfall=params%rainfall_available
  if (rainfall.gt.0.d0) geometry%precipitation(j)=geometry%z(i)/params%rainfall_height
  rainfall=max(0.d0,rainfall-geometry%precipitation(j))
  enddo
print*,params%rainfall_available,params%rainfall_height

geometry%precipitation=max(geometry%precipitation,params%rainfall_minimum)
stop

deallocate (network%donors,network%receiver,network%ndon)

call time_out ('find_precipitation')

return
end
