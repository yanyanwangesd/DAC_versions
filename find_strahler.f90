subroutine find_strahler (geometry,network,stack)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  integer i,j,k,h,l,maxs,occurance
  

  call time_in ('find_strahler')

  geometry%strahler=0
  do i=stack%nnode,1,-1
     j=stack%order(i)
     if (network%ndon(j).eq.0) then
        geometry%strahler(j) = 1
     else
        maxs = 0
        occurance = 0
        k = network%ndon(j)
        do h = 1,k
           l = network%donors(h,j)
           if (geometry%strahler(l).gt.maxs) then
              maxs = geometry%strahler(l)
              occurance = 1
           elseif (geometry%strahler(l).eq.maxs) then
              occurance = occurance +1
           endif
        enddo
        if (occurance.ge.2) then
           geometry%strahler(j) = maxs+1
        elseif(occurance.eq.1) then
           geometry%strahler(j) = maxs
        else
           print*, 'problem calculating strahler number'
        endif
     endif
  enddo
         
  call time_out ('find_strahler')

  return

end subroutine find_strahler
