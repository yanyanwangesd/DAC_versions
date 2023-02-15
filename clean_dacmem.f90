subroutine clean_dacmem (geometry,params,delaunay,network,stack,landscape)

  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay
  type (netw) network
  type (stck) stack
  type (parm) params
  type (lands) landscape


  deallocate (geometry%x,geometry%y,geometry%z)
  deallocate (geometry%xdiv,geometry%ydiv,geometry%zdiv)
  deallocate (geometry%u,geometry%v,geometry%w)
  deallocate (geometry%fix,geometry%surface, geometry%erosion_rate)
  deallocate (geometry%discharge,geometry%precipitation)
  deallocate (geometry%k,geometry%catchment)
  deallocate (geometry%nb)
  deallocate (geometry%nt)
  deallocate (geometry%nn)
  deallocate (geometry%nnodetri)
  deallocate (geometry%nndivtri)
  deallocate (geometry%nndivnode)
  deallocate (geometry%sediment_flux, geometry%chi)
  deallocate (geometry%strahler)
  deallocate (geometry%surface_share)
  if(params%transient_divide) then
     deallocate(geometry%erosion_rate_history)
  endif


   deallocate(delaunay%icon,delaunay%neighbours,delaunay%numdivides,delaunay%centers)

   deallocate(network%receiver,network%ndon,network%lakes,network%lakes_catch,network%donors)

   deallocate(stack%order)
return

if(params%f_varies_with_xyz)then
   deallocate (params%f_depends_on, params%f_variable_determined,params%f_polyseg, params%f_superpose  )  ! dimension f_num_sets
   deallocate ( params%f_timebound)       ! limit validity of condition block in time
   deallocate (params%poly )       ! polynomial coefficients
   deallocate ( params%pn )  ! degree/elements per line
endif

if(params%show_vtkfine)deallocate (landscape%xx,landscape%yy,landscape%zz,landscape%edot,landscape%fineicon)

end subroutine clean_dacmem
