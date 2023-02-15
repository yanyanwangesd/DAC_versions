subroutine find_discharge (geometry,network,stack)

  ! Routine to find the discharge (precipitation integrated over the catchment area) at each node
  ! Note that we also integrate erosion rate to obtain sediment flux

  ! To perform this operation one simply needs to go through the nodes in the stack in reverse order 

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  integer i,j,k

  call time_in ('find_discharge')

  geometry%discharge=geometry%surface*geometry%precipitation
  geometry%sediment_flux=geometry%surface*geometry%erosion_rate
  do i=stack%nnode,1,-1
     j=stack%order(i)
     k=network%receiver(j)
     if (k.ne.0) then
        geometry%discharge(k)=geometry%discharge(k)+geometry%discharge(j)
        geometry%sediment_flux(k)=geometry%sediment_flux(k)+geometry%sediment_flux(j)
     endif
  enddo

  call time_out ('find_discharge')

  return

end subroutine find_discharge
