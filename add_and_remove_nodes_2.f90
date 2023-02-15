subroutine add_remove_nodes(geometry,network,params,delaunay)

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay

  double precision l,xx,xd,fact3,fact2,lidon,lirec
  integer i,j,k,h 
  integer nrem_small_area,nrem_close_to_boundary,nrem_short_channel
  integer nadd_on_channel,nadd_between_channel,nadd_river_no_connection,nadd_boundary
  integer nnold, icount,flag,don,rec
  double precision xmin,xmax,ymin,ymax,min_dist
  integer num_samples, num_full_bin,samle_in_last_bin 
  double precision tau,ave_erosion 

  call time_in ('add_remove_nodes')


  !when advecting nodes there are special add and remove operations for which the following is needed:
  if (params%move_points) then 
     xmin = minval(geometry%x(1:geometry%nnode))
     xmax = maxval(geometry%x(1:geometry%nnode))
     ymin = minval(geometry%y(1:geometry%nnode))
     ymax = maxval(geometry%y(1:geometry%nnode))
  endif


  nnold=geometry%nnode !original number of nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! First stage: take care of all the cases where a new node should added !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! This section is to add nodes on channels where the channel connect two nodes that are not neighbors.
  nadd_river_no_connection=0

  do i=1,nnold ! loop over the nodes
     if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then
        flag=0
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
           if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
        enddo
        if(flag.eq.0)then !reciever is not a neighbor
           j=network%receiver(i)
           nadd_river_no_connection=nadd_river_no_connection+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0
           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
           ! print*, 'adding node', geometry%nnode, 'to prevent channel between non-neighboring node in:'
           ! print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
           ! Assign physical properties to new node as average of its endmembers
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0d0
           geometry%discharge(geometry%nnode)=0.0d0
           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
           geometry%sediment_flux(geometry%nnode)=0.0d0
           if (params%transient_divide) then
              do h=1,params%num_bins
                 geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                      geometry%erosion_rate_history(j,h))/2.d0
              enddo
           endif
        endif
     endif
  enddo
  !if (nadd_river_no_connection.gt.0) print*,  'nadd_river_no_connection', nadd_river_no_connection

  ! This section is to add nodes on channels. It checks whether nodes on the drainage network
  ! are too far apart or not. I adds nodes when the distance between two nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90) - SDW

  nadd_on_channel=0
  do i=1,nnold ! loop over the nodes
     if (network%receiver(i).ne.0) then ! only for connections that are part of the drainage network
        j=network%receiver(i)
        l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
        if (l.gt.params%lmax) then ! add a node on the connection
           nadd_on_channel=nadd_on_channel+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           ! Add a random perturbation in location of new node equal to 50% of the max length between nodes
           call random_number(xx)
           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                +(2.d0*xx-1.d0)*0.05d0*params%lmax 
           call random_number(xx)
           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                +(2.d0*xx-1.d0)*0.05d0*params%lmax
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
           !print*, 'adding node', geometry%nnode, 'in a too long channel in:'
           !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
           !note every physical properties will have to be updated at each new
           !node....
           ! Assign physical properties to new node as average of its endmembers
           ! added 15.6.2010 sdw  Not sure if these are all the node-local properties
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0d0
           geometry%discharge(geometry%nnode)=0.0d0
           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0           
           geometry%sediment_flux(geometry%nnode)=0.0d0 
           if (params%transient_divide) then
              do h=1,params%num_bins
                 geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                      geometry%erosion_rate_history(j,h))/2.d0
              enddo
           endif
        endif
     endif
  enddo


  !if (nadd_on_channel.gt.0) print*, 'nadd_on_channel',nadd_on_channel


  ! add node on non-channel divides that are longer than a given distance
  ! add node to lower existing node
  ! note that this loop does not take into acount boundary nodes.

  nadd_between_channel = 0
  ! the if statement below is needed to stop potentially crossing rivers -SDW
  ! (LG - maybe because the neighbor list was not updated)
  if(nadd_river_no_connection.eq.0.and.nadd_on_channel.eq.0)then
     do i=1,nnold ! loop over the original nodes
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
           j=geometry%nn(k,i)
           if (geometry%fix(i).eq.0.or.geometry%fix(j).eq.0) then !when both are boundary nodes treated with nadd_boundary
              if(network%receiver(i).ne.j.and.network%receiver(j).ne.i) then ! check non-channel
                 if(geometry%surface_share(k,i).gt.0.d0)then  ! only continue if divide is closer to j than i
                    l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
                    if (l.gt.params%ldivmax) then ! add node on the connection
                       xd=(geometry%surface_share(k,i)+1)*l/2.d0
                       if(xd.gt.params%xc)then  !  continue only if there is a fluvial segment to connection

                          nadd_between_channel=nadd_between_channel+1

                          geometry%nnode=geometry%nnode+1
                          !print*, 'adding node', geometry%nnode,'between',i,j
                          if (geometry%nnode.gt.geometry%nnode_max) then
                             print*,'problem while adding nodes'
                             STOP 'too many nodes added. Increase nnode_max'
                          endif
                          geometry%x(geometry%nnode)=geometry%x(i)+(xd-params%xc)/l*(geometry%x(j)-geometry%x(i))  
                          geometry%y(geometry%nnode)=geometry%y(i)+(xd-params%xc)/l*(geometry%y(j)-geometry%y(i)) 
                          if (params%transient_divide) then
                             if (params%h*params%m.eq.1) then
                                tau=log(l/params%xc)/ &
                                     (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))
                             else
                                tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                                     (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))/(1.d0-params%h*params%m)
                             endif
                             num_samples=floor(min(tau,params%time)/params%deltat)
                             if (num_samples.eq.0) then
                                fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                     (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                             else
                                if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                                   num_full_bin=params%num_bins
                                   samle_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                                else
                                   num_full_bin = INT(num_samples/params%sample_per_bin)+1
                                   samle_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                                endif
                                ave_erosion=0.d0
                                if (num_full_bin.eq.1) then
                                   ave_erosion=geometry%erosion_rate_history(i,1)
                                else
                                   ave_erosion=0.d0
                                   do h=1,num_full_bin-1
                                      ave_erosion=ave_erosion +geometry%erosion_rate_history(i,h)*params%sample_per_bin
                                   enddo
                                   ave_erosion=ave_erosion +geometry%erosion_rate_history(i,num_full_bin)*samle_in_last_bin
                                   ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + samle_in_last_bin)
                                endif
                                fact3=(ave_erosion/geometry%k(i)/ &
                                     (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                             endif
                          else
                             fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                          endif
                          if (params%hmn.eq.0.d0) then
                             fact2=log(params%xc/((xd)))                          
                             geometry%z(geometry%nnode)=geometry%z(i)-fact2*fact3 
                             if (fact3.eq.0.d0) then
                                call random_number(xx)           
                                geometry%z(geometry%nnode) = geometry%z(i) +xx*0.05d0
                                !print*, 'adding new node with not too different height',geometry%nnode,geometry%z(i),geometry%z(geometry%nnode)
                             endif
                          else
                             fact2=(1.d0/params%hmn)*(xd**(params%hmn)-(params%xc)**(params%hmn))
                             geometry%z(geometry%nnode)=geometry%z(i)+fact2*fact3 
                             if (fact3.eq.0.d0) then
                                call random_number(xx)           
                                geometry%z(geometry%nnode) = geometry%z(i) +xx*0.05d0
                                !print*, 'adding new node with the close height height',geometry%nnode,geometry%z(i),geometry%z(geometry%nnode)
                             endif
                          endif
                          !print*, 'adding node', geometry%nnode, 'because distance between channels is too large in:'
                          !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)

                          ! Assign physical properties to new node from neighbor node i
                          !  sdw  Not sure if these are all the node-local properties
                          geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                          geometry%surface(geometry%nnode)=params%amin
                          geometry%u(geometry%nnode)=geometry%u(i)
                          geometry%v(geometry%nnode)=geometry%v(i)
                          geometry%w(geometry%nnode)=geometry%w(i)
                          geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                          geometry%fix(geometry%nnode)=0
                          geometry%discharge(geometry%nnode)=0.0d0
                          geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                          geometry%k(geometry%nnode)=geometry%k(i)
                          geometry%sediment_flux(geometry%nnode)=0.0d0
                          network%ndon(i)=network%ndon(i)+1
                          network%ndon(geometry%nnode)=0
                          network%receiver(geometry%nnode)=i
                          if (params%transient_divide) then
                             do h=1,params%num_bins
                                geometry%erosion_rate_history(geometry%nnode,h)=geometry%erosion_rate_history(i,h)
                             enddo
                          endif
                       endif
                    endif
                 endif
              endif
           endif
        enddo
     enddo
  endif


  !if (nadd_between_channel.gt.0) print*, 'nadd_between_channel',nadd_between_channel

  ! This section is to add boundary nodes when the boundaries are advected. 
  ! It adds nodes when the distance between two boundary nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90)

  ! for some reason the nb array for boundary nodes is not trustable. Need to generate a new list.
  ! Each node will point to its neighbor in anticlockwise direction.

  !divided into four section according tobottom, right, top and left boundaries

  ! try again with the neighbors list 1/4/11
  nadd_boundary=0
  if (params%move_points) then
     do i=1,nnold
        if (geometry%fix(i).ge.1) then
           if (geometry%y(i).eq.ymin) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%y(j).eq.ymin.and.geometry%x(j).gt.geometry%x(i)) then !this is the closesent neighbor to the right
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%y(geometry%nnode)=geometry%y(i)
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0 !should be zero
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large in:'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=geometry%fix(i)
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%x(i).eq.xmax) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%x(j).eq.xmax.and.geometry%y(j).gt.geometry%y(i)) then !this is the closesent neighbor to the top
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%x(geometry%nnode)=geometry%x(i)
                       call random_number(xx)              
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0

                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmax and yi=ymin)
                       !geometry%fix(geometry%nnode)=geometry%fix(i)
                       if (geometry%y(j) .eq. ymax) then
                           geometry%fix(geometry%nnode)=geometry%fix(i)
                       else
                           geometry%fix(geometry%nnode)=geometry%fix(j)
                       endif
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmax and yi=ymin)

                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%y(i).eq.ymax) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%y(j).eq.ymax.and.geometry%x(j).lt.geometry%x(i)) then !this is the closesent neighbor to the left
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%y(geometry%nnode)=geometry%y(i)
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0 !should be zero
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=geometry%fix(i)
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%x(i).eq.xmin) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%x(j).eq.xmin.and.geometry%y(j).lt.geometry%y(i)) then !this is the closesent neighbor to the bottom
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%x(geometry%nnode)=geometry%x(i)
                       call random_number(xx)              
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0

                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmin and yi=ymin)
                       !geometry%fix(geometry%nnode)=geometry%fix(i)
                       if (geometry%y(j) .eq. ymin) then
                           geometry%fix(geometry%nnode)=geometry%fix(i)
                       else
                           geometry%fix(geometry%nnode)=geometry%fix(j)
                       endif
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmin and yi=ymin)

                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
        endif
     enddo
  endif


 ! if (nadd_boundary.gt.0) print*, 'nadd_boundary',nadd_boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Second stage: take care of all the cases where a new node should removed !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! this section remove nodes based on minimum surface area in headwater
  nrem_small_area=0
  do i=1, nnold
     !because we remove nodes, need to check if we are still within the bounds of ald nodes. 
     if(i.le.geometry%nnode-nadd_on_channel-nadd_between_channel-nadd_river_no_connection-nadd_boundary)then 
        if(geometry%fix(i).eq.0.and.network%ndon(i).eq.0.and.geometry%surface(i).lt.params%amin)then
           nrem_small_area=nrem_small_area+1
           geometry%nnode=geometry%nnode-1  
           !print*, 'removing node', i, 'because its a leaf with too small drainage area'
           !print*, geometry%x(i), geometry%y(i), geometry%z(i)
           do j=1,i-1
              if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
           enddo
           if(i.le.geometry%nnode)then
              do j=i,geometry%nnode
                 geometry%x(j)=geometry%x(j+1)
                 geometry%y(j)=geometry%y(j+1)
                 geometry%z(j)=geometry%z(j+1)
                 geometry%u(j)=geometry%u(j+1)
                 geometry%v(j)=geometry%v(j+1)
                 geometry%w(j)=geometry%w(j+1)
                 geometry%fix(j)=geometry%fix(j+1)
                 geometry%surface(j)=geometry%surface(j+1)
                 geometry%discharge(j)=geometry%discharge(j+1)
                 geometry%precipitation(j)=geometry%precipitation(j+1)
                 geometry%k(j)=geometry%k(j+1)
                 geometry%nb(j)=geometry%nb(j+1)
                 geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                 geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                 network%ndon(j)=network%ndon(j+1)
                 if (network%receiver(j+1).gt.i) then
                    network%receiver(j)=network%receiver(j+1)-1
                 else
                    network%receiver(j)=network%receiver(j+1)
                 endif
              enddo
              if (params%transient_divide) then
                 do j=i,geometry%nnode
                    geometry%erosion_rate_history(j,:)=geometry%erosion_rate_history(j+1,:)
                 enddo
              endif
           endif
        endif
     endif
  enddo

!  if (nrem_small_area.gt.0) print*, 'nrem_small_area', nrem_small_area

!!! the next remove section operates also on nodes that were just now added. 
!!!For that reason the donor receiver array needs to be updated
  if (nrem_small_area.gt.0.or.nadd_on_channel.gt.0.or.nadd_between_channel.gt.0&
       &.or.nadd_river_no_connection.gt.0.or.nadd_boundary.gt.0)then
     network%nnode=geometry%nnode
     network%donors=0
     network%ndon=0
     do i=1,network%nnode
        k=network%receiver(i)
        if (k.ne.0) then
           network%ndon(k)=network%ndon(k)+1
           network%donors(network%ndon(k),k)=i
        endif
     enddo
  endif

  ! this section remove nodes that defines too short channel links mostly resulting from nnode-nadd_on_channel
  ! removes nodes with only one donor, such that the distance between its donor and receiver is smaller than 0.5*params%lmax 
  nrem_short_channel=0

 if (params%move_points.and.maxval(geometry%nb).ge.geometry%nnmax-1) then 
    do i=1,nnold
     !because we remove nodes, need to check if we are still within the bounds of ald nodes. 
553    if(i.le.geometry%nnode+nrem_small_area-nadd_on_channel-nadd_between_channel-nadd_river_no_connection-nadd_boundary)then        
          if(geometry%fix(i).eq.0.and.network%ndon(i).eq.1.and.network%receiver(i).ne.0.and.geometry%nb(i).eq.4) then ! the nb conditions insures that this is not a new node
             don=network%donors(1,i)
             rec=network%receiver(i)
             lidon=dsqrt((geometry%x(i)-geometry%x(don))**2+(geometry%y(i)-geometry%y(don))**2) 
             lirec=dsqrt((geometry%x(i)-geometry%x(rec))**2+(geometry%y(i)-geometry%y(rec))**2)
             if  (lidon+lirec.lt.0.75*params%lmax) then
                !print*, 'removing node', i, 'because its defining too short channel segments'
                nrem_short_channel=nrem_short_channel+1
                geometry%nnode=geometry%nnode-1            
                ! first reorganize network
                network%receiver(don) = rec ! now h is a donor of k
                do j=1,i-1
                   if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
                enddo
                if(i.le.geometry%nnode)then !not the last node
                   do j=i,geometry%nnode
                      geometry%x(j)=geometry%x(j+1)
                      geometry%y(j)=geometry%y(j+1)
                      geometry%z(j)=geometry%z(j+1)
                      geometry%u(j)=geometry%u(j+1)
                      geometry%v(j)=geometry%v(j+1)
                      geometry%w(j)=geometry%w(j+1)
                      geometry%fix(j)=geometry%fix(j+1)
                      geometry%surface(j)=geometry%surface(j+1)
                      geometry%discharge(j)=geometry%discharge(j+1)
                      geometry%precipitation(j)=geometry%precipitation(j+1)
                      geometry%k(j)=geometry%k(j+1)
                      geometry%nb(j)=geometry%nb(j+1)
                      geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                      geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                      network%ndon(j)=network%ndon(j+1)
                      if (network%receiver(j+1).gt.i) then
                         network%receiver(j)=network%receiver(j+1)-1
                      else
                         network%receiver(j)=network%receiver(j+1)
                      endif
                   enddo
                   if (params%transient_divide) then
                      do j=i,geometry%nnode
                         geometry%erosion_rate_history(j,:)=geometry%erosion_rate_history(j+1,:)
                      enddo
                   endif
                endif
                !need to continuously update donors and receivers
                network%nnode=geometry%nnode
                network%donors=0
                network%ndon=0
                do j=1,network%nnode
                   k=network%receiver(j)
                   if (k.ne.0) then
                      network%ndon(k)=network%ndon(k)+1
                      network%donors(network%ndon(k),k)=j
                   endif
                enddo
                go to 553 !without increasing i becuase the node array is smaller.
             endif
          endif
       endif
    enddo
 endif
 




             

  !if (nrem_short_channel.gt.0) print*, 'nrem_short_channel', nrem_short_channel

!!! the next remove section operates also on nodes that were just now added. 
!!!For that reason the donor receiver array needs to be updated
  if (nrem_short_channel.gt.0.or.nrem_small_area.gt.0.or.nadd_on_channel.gt.0&
       &.or.nadd_between_channel.gt.0.or.nadd_river_no_connection.gt.0.or.nadd_boundary.gt.0)then
     network%nnode=geometry%nnode
     network%donors=0
     network%ndon=0
     do i=1,network%nnode
        k=network%receiver(i)
        if (k.ne.0) then
           network%ndon(k)=network%ndon(k)+1
           network%donors(network%ndon(k),k)=i
        endif
     enddo
  endif


  ! This section removes nodes that are too close to the boundary due to advection
  ! Both internal modes and boundary nodes

  nrem_close_to_boundary = 0
  if (params%move_points) then 
     min_dist =  params%max_adv*params%deltat !if a distance of a node to the boundary < min_dist-->remove it
     do i=1, geometry%nnode !also new nodes that were just added
554     if (i.le.geometry%nnode) then ! i can grow above geometry%nnode because geometry%nnode can become smaller duringt this loop
           if(geometry%fix(i).eq.0.or.& !internal nodes
                &geometry%y(i).eq.ymin.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not coreners
                &geometry%y(i).eq.ymax.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not coreners
                &geometry%x(i).eq.xmin.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax.or.& !boundary nodes but not coreners
                &geometry%x(i).eq.xmax.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax) then!boundary nodes but not coreners
              if (abs(geometry%x(i)-xmin).lt.min_dist.and.geometry%x(i).ne.xmin.or.& !close to left
                   &abs(geometry%x(i)-xmax).lt.min_dist.and.geometry%x(i).ne.xmax.or.& !close to right
                   &abs(geometry%y(i)-ymin).lt.min_dist.and.geometry%y(i).ne.ymin.or.& !close to bottom
                   &abs(geometry%y(i)-ymax).lt.min_dist.and.geometry%y(i).ne.ymax) then !close to top
                 !print*, 'removing node', i, 'because its too close to a boundary'
                 !print*, geometry%x(i), geometry%y(i), geometry%z(i) 
                 !print*, 'updated number of nodes is',geometry%nnode
                 nrem_close_to_boundary = nrem_close_to_boundary +1
                 geometry%nnode=geometry%nnode-1                 
                 ! first reorganize network
                 k = network%receiver(i)
                 if (k.ne.0) then !i is not a lake
                    do j = 1,network%ndon(i)
                       h = network%donors(j,i) ! h is a donor of i
                       network%receiver(h) = k ! now h is a donor of k
                    enddo
                 else! i is a lake
                    do j = 1,network%ndon(i)
                       h = network%donors(j,i) ! h is a donor of i
                       network%receiver(h) = 0 ! now h is a lake
                    enddo
                 endif
                 do j=1,i-1
                    if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
                 enddo
                 if(i.le.geometry%nnode)then !not the last node
                    do j=i,geometry%nnode
                       geometry%x(j)=geometry%x(j+1)
                       geometry%y(j)=geometry%y(j+1)
                       geometry%z(j)=geometry%z(j+1)
                       geometry%u(j)=geometry%u(j+1)
                       geometry%v(j)=geometry%v(j+1)
                       geometry%w(j)=geometry%w(j+1)
                       geometry%fix(j)=geometry%fix(j+1)
                       geometry%surface(j)=geometry%surface(j+1)
                       geometry%discharge(j)=geometry%discharge(j+1)
                       geometry%precipitation(j)=geometry%precipitation(j+1)
                       geometry%k(j)=geometry%k(j+1)
                       geometry%nb(j)=geometry%nb(j+1)
                       geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                       geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                       network%ndon(j)=network%ndon(j+1)
                       if (network%receiver(j+1).gt.i) then
                          network%receiver(j)=network%receiver(j+1)-1
                       else
                          network%receiver(j)=network%receiver(j+1)
                       endif
                    enddo
                    if (params%transient_divide) then
                       do j=i,geometry%nnode
                          geometry%erosion_rate_history(j,:)=geometry%erosion_rate_history(j+1,:)
                       enddo
                    endif
                 endif
                 !need to continuously update donors and receivers
                 network%nnode=geometry%nnode
                 network%donors=0
                 network%ndon=0
                 do j=1,network%nnode
                    k=network%receiver(j)
                    if (k.ne.0) then
                       network%ndon(k)=network%ndon(k)+1
                       network%donors(network%ndon(k),k)=j
                    endif
                 enddo
                 go to 554 !without increasing i becuase the node array is smaller.
              endif
           endif
        endif
     enddo
  endif
 ! if (nrem_close_to_boundary.gt.0) print*,  'nrem_close_to_boundary',  nrem_close_to_boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Short debuging section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(network%nnode.ne.geometry%nnode)print*, 'error in addremovenodes network nnode wrong'

  icount=0
  do i=1,geometry%nnode
     if(geometry%fix(i).ge.1)then
        icount=icount+1
        go to 555
     endif
     if(network%receiver(i).eq.0)then
        icount=icount+1
        go to 555
     endif
     j=network%receiver(i)
     if(j.gt.0.and.j.le.geometry%nnode)then
        icount=icount+1
        go to 555
     endif
     print*,' node not in network: ',i, network%receiver(i)
555  continue
  enddo



  call time_out ('add_remove_nodes')

  return 
end subroutine add_remove_nodes
