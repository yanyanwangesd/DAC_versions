subroutine calculate_delaunay (geometry,delaunay)

  ! interface routine to the delaun routine written by Malcolm Sambridge
  ! calculates the Delaunay triangulation around a set of points
  ! the number of triangles is stored in delaunay%ntriangles
  ! the triangulation is stored in delaunay%icon and the 
  ! neighboouring triangle information in delaunay%neighbours

  ! This is a "utility" routine and should not be changed

  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay

  real*8,dimension(:,:),allocatable::points
  real*8,dimension(:,:),allocatable::centers
  integer,dimension(:,:),allocatable::e,v
  integer,dimension(:),allocatable::vis_tlist,vis_elist,add_tlist
  integer num,numtri_max,numtri,nv_max,del_mode,nfirst,itstart
  logical,dimension(:),allocatable::inactive,subset
  real*8 eps
  integer i, j
  real*8 xy(2,3), ctr(2)
  integer node

  call time_in ('calculate_delaunay')

  num=geometry%nnode
  numtri_max=num*3
  nv_max=num
  del_mode=0
  nfirst=1
  itstart=1
  eps=1.d-10

  allocate (points(2,num))
  allocate (centers(3,numtri_max))
  allocate (e(3,numtri_max),v(3,numtri_max))
  allocate (vis_tlist(nv_max),vis_elist(nv_max),add_tlist(nv_max))
  allocate (inactive(num),subset(num))


  do i=1,geometry%nnode
     points(1,i)=geometry%x(i)
     points(2,i)=geometry%y(i)
  enddo
  inactive=.FALSE.
  subset=.TRUE.

  !print*,'going to call delaun'
  call delaun (points,num,e,v,numtri,numtri_max, &
       vis_tlist,vis_elist,add_tlist,eps,nv_max, &
       del_mode,inactive,nfirst,itstart,subset)
  !print*,'after delaun'

  if (delaunay%ntriangles.ne.0) deallocate (delaunay%icon,delaunay%neighbours,delaunay%centers, delaunay%numdivides)

  delaunay%ntriangles=numtri
  allocate (delaunay%icon(3,delaunay%ntriangles),delaunay%neighbours(3,delaunay%ntriangles))
  allocate (delaunay%numdivides(4,delaunay%ntriangles),delaunay%centers(3,delaunay%ntriangles))
  delaunay%centers=0
  delaunay%numdivides=0

  delaunay%icon=v(:,1:numtri)
  delaunay%neighbours=e(:,1:numtri)

  ! call ccentres(points,v,numtri,centers)  ! calculates centres of all Delaunay circumcircles
  ! delaunay%centers(1:2,1:numtri) = centers(1:2,1:numtri)
  ! To calculate the centroids of the Delaunay

  do i = 1,delaunay%ntriangles
     xy(1,1) = geometry%x(delaunay%icon(1,i))
     xy(2,1) = geometry%y(delaunay%icon(1,i))
     xy(1,2) = geometry%x(delaunay%icon(2,i))
     xy(2,2) = geometry%y(delaunay%icon(2,i))
     xy(1,3) = geometry%x(delaunay%icon(3,i))
     xy(2,3) = geometry%y(delaunay%icon(3,i))
     call calculate_centroid (xy, ctr)
     delaunay%centers(1:2,i) = ctr
     !   print*, 'center of mass is', ctr(1), ctr(2)
     delaunay%centers(3,i) = 0.
     delaunay%numdivides(1,i) = 0
  enddo

  geometry%nnodetri=0
  geometry%nt=0
  do i= 1,delaunay%ntriangles
     do j=1,3
        node=delaunay%icon(j,i)
        geometry%nt(node)=geometry%nt(node)+1
        geometry%nnodetri(geometry%nt(node),node)=i
     enddo
  enddo


  deallocate (points,centers,e,v,vis_tlist,vis_elist,add_tlist,inactive,subset)

  call time_out ('calculate_delaunay')

  return

end subroutine calculate_delaunay
