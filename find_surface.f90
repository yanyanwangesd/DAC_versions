subroutine find_surface (geometry,network,params)

! Utility routine to compute the surface area associated with each node
! Uses Lasserre's recurrence algorithm
! Should not be changed

use definitions

implicit none

type (geom) geometry
type (netw) network
type (parm) params
integer n,i,j,k,m
real,dimension(:,:),allocatable::pp,a
real,dimension(:),allocatable::xx,b
double precision xmin,xmax,ymin,ymax,surfscale,dsurf
double precision,dimension(:),allocatable::dsurface
real vol

call time_in ('find_surface')

allocate (xx(2),pp(2,geometry%nnmax),a(geometry%nnmax,2),b(geometry%nnmax))

xmin=minval(geometry%x)
xmax=maxval(geometry%x)
ymin=minval(geometry%y)
ymax=maxval(geometry%y)
surfscale=(xmax-xmin)*(ymax-ymin)/geometry%nnode

n=2
  do i=1,geometry%nnode
  xx(1)=geometry%x(i)
  xx(2)=geometry%y(i)
    do j=1,geometry%nb(i)
    pp(1,j)=geometry%x(geometry%nn(j,i))
    pp(2,j)=geometry%y(geometry%nn(j,i))
    enddo
  m=geometry%nb(i)
  call first_voronoi (xx,pp,n,m,geometry%nnmax,2,a,b,vol)
    if (vol.gt.surfscale*10.) then
    geometry%surface(i)=surfscale
    elseif (vol.gt.0.) then
    geometry%surface(i)=vol
    else
    geometry%surface(i)=surfscale
    endif
  enddo

deallocate (xx,pp,a,b)

! adjust surface from divide calculations
! one assumes that each node has the potential to give an equal share of its surface
! to each of its neighbours; how much of that share is given depends on the
! value of the surface_share parameter.
! if surface_share is 0 nothing is shared
! if surface_share is > 0 node i takes a share of the share of j proportional to surface_share
! if surface_share is < 0 node i gives a share of its share proportional to -surface_share

if (params%divide.or.params%small_divide) then
allocate (dsurface(geometry%nnode))
dsurface=0.d0
  do i=1,geometry%nnode
    do k=1,geometry%nb(i)
    j=geometry%nn(k,i)
      if (geometry%surface_share(k,i).gt.0.d0) then
      dsurf=geometry%surface(j)/geometry%nb(j)*geometry%surface_share(k,i)
      dsurface(i)=dsurface(i)+dsurf
      dsurface(j)=dsurface(j)-dsurf
      elseif (geometry%surface_share(k,i).lt.0.d0) then
      dsurf=geometry%surface(i)/geometry%nb(i)*geometry%surface_share(k,i)
      dsurface(i)=dsurface(i)+dsurf
      dsurface(j)=dsurface(j)-dsurf
      endif
    enddo
  enddo
 
  do i=1,geometry%nnode
   geometry%surface(i)=geometry%surface(i)+dsurface(i)
  enddo 
deallocate (dsurface)
endif

call time_out ('find_surface')

return

end subroutine find_surface
