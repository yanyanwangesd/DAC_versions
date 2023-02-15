subroutine find_catchment (geometry,network,stack,delaunay,params)

  ! Routine to compute catchment number (name) for each node
  ! This name (number) is taken as the number of the baselevel node
  ! defining the catchment

  ! To perform this operation, one simply needs to go through the stack resetting the
  ! catchemnt name each time a no receiver node is found

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (del) delaunay
  type (parm) params

  integer i,j,jj,icatch,k,jk,kk,ilcatch
  integer ilk,np,ispill,jkd
  double precision zmin,xx

  call time_in ('find_catchment')

  do i=1,stack%nnode
     j=stack%order(i)
     if (network%receiver(j).eq.0) then
        icatch=j
        ilcatch=1
        if (geometry%fix(j).ge.1) ilcatch=0
     endif
     geometry%catchment(j)=icatch
     network%lakes_catch(j)=ilcatch
  enddo

  call time_out ('find_catchment')
  return

  ! Below is the exclude local lakes 
  ! It first finds lowest point on the edge of the catchment and use this as a
  ! location for a sill (or spillway..., depending on your own version of English
  ! ;))

  ! this does not really work (it should be improved/tested)

  if (network%nlake.gt.0) then

     do i=1,network%nlake
        zmin=maxval(geometry%z)
        ispill=stack%order(network%lakes(i))
        j=network%lakes(i)
        icatch=geometry%catchment(stack%order(j))
        jk=stack%order(j+1)
        do while (geometry%catchment(jk).eq.icatch)
           do kk=1,geometry%nb(jk)
              if (geometry%catchment(geometry%nn(kk,jk)) .ne. icatch ) then
                 if (geometry%z(geometry%nn(kk,jk)) .lt. geometry%z(jk)) then
                    if (geometry%z(jk) .lt. zmin) then
                       if (network%lakes_catch(geometry%nn(kk,jk)).lt.1) then
                          zmin=geometry%z(jk)
                          ispill=geometry%nn(kk,jk)
                          jkd=jk
                       endif
                    endif
                 endif
              endif
           enddo
           j=j+1
           jk=stack%order(j)
        enddo

        if (icatch.ne.geometry%catchment(ispill)) then
           call random_number (xx)
           geometry%z(jkd)=geometry%z(ispill)+xx/100.
           network%receiver(stack%order(network%lakes(i)))=jkd
           network%receiver(jkd)=ispill
           !       network%ndon(ispill)=network%ndon(ispill)+1
           !       if (network%ndon(ispill).gt.network%ndonmax) stop 'network%ndonmax too small, in find_catchments'
           !       network%donors(network%ndon(ispill),ispill)=stack%order(network%lakes(i))
        endif

     enddo

     if (network%nlake.ne.0) then
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

     call find_order (geometry,network,stack,delaunay,params) 

     do k=1,stack%nnode
        j=stack%order(k)
        if (network%receiver(j).eq.0) icatch=j
        geometry%catchment(j)=icatch
     enddo

  endif


  call time_out ('find_catchment')

  return

end subroutine find_catchment
