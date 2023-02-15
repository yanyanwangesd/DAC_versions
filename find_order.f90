subroutine find_order (geometry,network,stack,delaunay,params)

  ! ROutine to compute the order in which the nodes should be processed
  ! to compute the erosion and (in reverse) compute the discharge

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (del) delaunay
  type (parm) params

  integer i,j,nnodecheck,ii,k
  double precision zmn,xx

  call time_in ('find_order')

  if (stack%nnode.eq.0.or.params%add_nodes) then
     allocate (stack%order(geometry%nnode))
  endif

  ! goes through the list of nodes and starts building the stack for each node that has no receiver
  ! ie, it is a local minimum or a base level (fix=1) node
  ! also check whether local minimum is in the middle of the grid (i.e. local
  ! lake)

  stack%nnode=0

  network%nlake=0
  network%lakes=0

  do i=1,geometry%nnode
     if (network%receiver(i).eq.0.or.geometry%fix(i).ge.1) then
        !      if (geometry%fix(i).eq.0) then
        !        network%nlake=network%nlake+1
        !        network%lakes(network%nlake)=stack%nnode+1
        !zmn=geometry%z(i)*1000.
        !do k=geometry%nb(i),1,-1 ! loop over all neighbours
        !  j=geometry%nn(k,i)
        !  if (geometry%z(j).lt.zmn) zmn=geometry%z(j)
        !enddo
        !xx=rand()
        !geometry%z(j)=zmn !+xx
        !      endif
        nnodecheck=stack%nnode
        call add_to_stack (i,stack,network)
        if (nnodecheck.eq.stack%nnode) then
           print*,i,'this node is not in the stack'
           print*,geometry%x(i),geometry%y(i),geometry%z(i),network%receiver(i),network%ndon(i)
           pause
        endif
     endif
  enddo

  !print*,network%nlake,'lakes'
  if (stack%nnode.ne.geometry%nnode) then
     print*,stack%nnode,geometry%nnode,'stack nnode and geometry nnode are different...why?'
     print*,'TOPO min-max',minval(geometry%z),maxval(geometry%z)
     !call show (geometry,delaunay,network,stack,params)
     !pause
     do i=1, geometry%nnode
        do j=1, stack%nnode
           if(i.eq.stack%order(j))go to 111
        enddo
        print*,'missing node in stack', i
111     continue
     enddo
  endif
  !print*,stack%nnode,geometry%nnode

  call time_out ('find_order')

  return

end subroutine find_order

! ---

recursive subroutine add_to_stack (i,stack,network)

! Recursive routine to go through the node upstream and keeping track of confluences

use definitions

implicit none

type (stck) stack
type (netw) network
integer i,j

stack%nnode=stack%nnode+1
if (stack%nnode.gt.size(stack%order)) then
print*,stack%nnode
stop 'error in add_to_stack'
endif
stack%order(stack%nnode)=i

if (network%ndon(i).eq.0) return

  do j=1,network%ndon(i)
  call add_to_stack (network%donors(j,i),stack,network)
  enddo

return

end subroutine add_to_stack
