subroutine find_network (geometry,delaunay,network,iflag)

  ! Routine to extract from the triangulation the list of neighbour nodes,
  ! donor nodes and receiver node for each node

  ! This is a utility routine and should not be changed

  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay
  type (netw) network
  integer nnode,it,i,j,k,i1,i2,iflag, bflag
  double precision slope,slopemax

  call time_in ('create_network')

  if (iflag.eq.0) then
     if (network%nnode.ne.0) deallocate (network%donors,network%receiver,&
          network%lakes_catch,network%lakes,network%ndon)

     network%nnode=geometry%nnode
     network%ndonmax=geometry%nnmax
     allocate (network%donors(network%ndonmax,geometry%nnode_max))
     allocate (network%receiver(geometry%nnode_max))
     allocate (network%lakes_catch(geometry%nnode_max))
     allocate (network%lakes(geometry%nnode_max))
     allocate (network%ndon(geometry%nnode_max))

  endif

  ! Go through triangles and compute neighbour list

  nnode=network%nnode
  geometry%nb=0
  geometry%nn=0
  do it=1,delaunay%ntriangles
     do k=1,3
        i1=delaunay%icon(k,it)
        i2=delaunay%icon(mod(k,3)+1,it)
        geometry%nb(i1)=geometry%nb(i1)+1
        if (geometry%nb(i1).gt.geometry%nnmax)  stop 'geometry%nnmax must be increased'
        geometry%nn(geometry%nb(i1),i1)=i2
        if (delaunay%neighbours(mod(k+1,3)+1,it).eq.0) then
           geometry%nb(i2)=geometry%nb(i2)+1
           if (geometry%nb(i2).gt.geometry%nnmax) stop 'geometry%nnmax must be increased'
           geometry%nn(geometry%nb(i2),i2)=i1
        endif
     enddo
  enddo

  

  if (iflag.eq.1) return

  ! Using the steepest slope among the neighbouring nodes
  ! compute the receiver node of each node (unique)

  network%receiver=0
  do i=1,geometry%nnode
     if (geometry%nb(i).eq.0) stop 'Disconnected node in create_network'
     slopemax=0.d0
     do j=1,geometry%nb(i)
        k=geometry%nn(j,i)
        slope=(geometry%z(i)-geometry%z(k))
        bflag=0
        if(geometry%fix(k).ne.0)  bflag = geometry%boundary( geometry%fix(k) )
        if (slope.gt.0.d0) then
           slope=slope/sqrt((geometry%x(i)-geometry%x(k))**2+(geometry%y(i)-geometry%y(k))**2)
           if (slope.gt.slopemax     .and. bflag.ge.0) then
              slopemax=slope
              network%receiver(i)=k
           endif
        endif
     enddo
  enddo

  ! Using the receiver node for each node, compute the donors nodes for each node

  network%donors=0
  network%ndon=0
  do i=1,network%nnode
     k=network%receiver(i)
     if (k.ne.0) then
        network%ndon(k)=network%ndon(k)+1
        if (network%ndon(k).gt.network%ndonmax) stop 'network%ndonmax too small'
        network%donors(network%ndon(k),k)=i
     endif
  enddo

  call time_out ('create_network')

  return

end subroutine find_network
