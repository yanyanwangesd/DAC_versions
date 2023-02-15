subroutine update_erosion_rate_history (geometry,params)

  ! subroutine to find the average erosion rate of each node for the 
  ! calculation of the transient elevation of divide



  use definitions

  implicit none

  type (geom) geometry
  type (parm) params
  integer i,j,num_full_bin, sample_in_last_bin


  call time_in ('update_erosion_rate_history')

  if ((params%num_bins-1)*params%sample_per_bin.lt.params%istep) then
     num_full_bin=params%num_bins
     sample_in_last_bin=params%istep-(params%num_bins-1)*params%sample_per_bin
  else
     num_full_bin = INT(params%istep/params%sample_per_bin)+1
     sample_in_last_bin=params%istep-(num_full_bin-1)*params%sample_per_bin
  endif
  do i=1,geometry%nnode
     if (num_full_bin.ne.1) then
        !first treat the last full bin
        geometry%erosion_rate_history(i,num_full_bin)= &
             (geometry%erosion_rate_history(i,num_full_bin)*sample_in_last_bin + &
             geometry%erosion_rate_history(i,num_full_bin-1))/(sample_in_last_bin+1)
        !then treat the rest of the bins
        do j=num_full_bin-1,2,-1
           geometry%erosion_rate_history(i,j)= &
                (geometry%erosion_rate_history(i,j)*(params%sample_per_bin-1) + &
                geometry%erosion_rate_history(i,j-1))/(params%sample_per_bin)
        enddo
        !last treat the first bin
        geometry%erosion_rate_history(i,1)= &
             (geometry%erosion_rate_history(i,1)*(params%sample_per_bin-1) + &
             geometry%erosion_rate(i))/(params%sample_per_bin)
     else
        geometry%erosion_rate_history(i,1)= &
             (geometry%erosion_rate_history(i,1)*(params%istep-1) + &
             geometry%erosion_rate(i))/(params%istep)
     endif

!!$     if (i.eq.41) then
!!$        print*, 'info about',i
!!$        do j=1,num_full_bin
!!$           print*, params%istep
!!$           print*, num_full_bin,sample_in_last_bin 
!!$           print*, 'geometry%erosion_rate_history(i)',geometry%erosion_rate_history(i,j)
!!$        enddo
!!$
!!$     endif


  enddo




  call time_out ('update_erosion_rate_history')

  return

end subroutine update_erosion_rate_history
