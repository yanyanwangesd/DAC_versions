program DivideAndCapture

  ! Beta version of new cascade algorithm that incorporates river capture
  ! but not the computation of precise divide locations
  ! Jean BRAUN 20 july 2009

  use definitions

  implicit none

  ! for a definition of the derived types see the file module_definitions.f90

  type (geom) geometry
  type (del) delaunay
  !type (del) rddelaunay
  type (netw) network
  type (stck) stack
  type (timr) timer
  type (parm) params
  type (lands)  landscape	
  character*11  title
  integer old_ncapture, i, j, k, m, l, n
  double precision li,lj
  ! sets a common for timing of various subroutines
  ! the results are printed to the screen at the end of the run

  common /global_timer/ timer

  call start_timing
  call floader( params, geometry)
  call initialize_parameters (params) ! initializes the model parameters
  if (params%read_restart) then 
     call reinit_geometry(geometry,params,network,delaunay)
  else
     call initialize_geometry (geometry,params) ! initializes the model geometry and history
  endif

  call set_params_as_function ( geometry, params )
 
  if (.not.params%read_restart) then
     delaunay%ntriangles=0
     call calculate_delaunay (geometry,delaunay) ! calculates Delaunay and neighbour (triangle) information
  endif
 
  if (params%read_restart) then
     call find_network (geometry,delaunay,network,1) ! find neighbour list
  else
     network%nnode=0
     call find_network (geometry,delaunay,network,0) ! find neighbour list and donor/receiver list
  endif
 
  if (.not.params%read_restart) then
     call find_precipitation (geometry,params) ! orography (not fully working yet)
     call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point
  endif

  call reboundary(geometry, params)

  stack%nnode=0
  call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms

  call find_strahler(geometry,network,stack)
 
  call output_z_tau(params,geometry,network,stack)

  call VTK (geometry,delaunay,params,network,stack)

  if (params%show_vtkfine)  call finetri (geometry, params,delaunay, network, landscape)
  if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)

  if (params%ascii) call write_ascii(geometry,params,network)
  call output_geometry(geometry,params,network,delaunay)


  !more frequent global data on the system is output to these files
  if (params%ascii) then
      open (80,file='ASCII/capture_data',status='unknown')
      write(80,'(F8.1,i10)') params%time, geometry%ncapture
      old_ncapture = geometry%ncapture
  endif

  do while (params%time.lt.params%tfinal) ! start of time loop
     
    
    
     params%deltat=min(params%deltat,params%tfinal-params%time)
     params%time=params%time+params%deltat
     params%istep=params%istep+1

     call set_params_as_function ( geometry, params )


     call uplift_and_advect (geometry,params) ! move and uplift the grid according to the prescribed velocity field


     if (params%move_points) then !with horizontal advection
        call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation        
        call find_network (geometry,delaunay,network,1) ! find neighbour list only       
     endif

     if (params%add_nodes) then 
        call add_remove_nodes(geometry,network,params,delaunay)       
        call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation      
        call find_network (geometry,delaunay,network,1) ! find neighbour list only         
     endif
     if (params%transient_divide) call update_erosion_rate_history(geometry,params)

     call captures_and_divides (geometry,network,params,delaunay) ! find potential captures and adjusts network accordingly  
    
     call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point
     if (params%add_nodes) deallocate(stack%order)
     call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms
     call find_catchment (geometry,network,stack,delaunay,params) ! find catchments
     call find_discharge (geometry,network,stack) ! compute discharge (improved cascade algorithm)
     call erode (geometry,network,stack,delaunay,params) ! erode
     
     if ((params%istep/params%freq)*params%freq.eq.params%istep) then	
        call find_strahler(geometry,network,stack) 
         if (params%show_vtkfine) call finetri (geometry, params,delaunay,network, landscape)
        call output_z_tau(params,geometry,network,stack)
        call VTK (geometry,delaunay,params,network,stack)
        if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)
        if (params%ascii) call write_ascii(geometry,params,network)
        call output_geometry(geometry,params,network,delaunay)
        print*, 'max number of neighbors is:', maxval(geometry%nb)
     endif

     !frequent data is output here
     if (params%ascii) then
         if (old_ncapture.ne.geometry%ncapture.or.MOD(params%istep,100).eq.0) then
            write(80,'(e15.7,i10)') params%time, geometry%ncapture
            old_ncapture = geometry%ncapture
         endif
     endif

  enddo ! end of time loop

  if (params%ascii) then
    write(80,'(e15.7,i10)') params%time, geometry%ncapture
    close (80)
  endif

  call find_strahler(geometry,network,stack)
  if (params%show_vtkfine)  call finetri (geometry, params,delaunay,network, landscape)
  call output_z_tau(params,geometry,network,stack)
  call VTK (geometry,delaunay,params,network,stack)
  if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)
  if (params%ascii) call write_ascii(geometry,params,network)
  call output_geometry(geometry,params,network,delaunay)
  call find_hack(geometry,network,stack)
  call show_timing

  call clean_dacmem(geometry, params,delaunay,network,stack, landscape)

end program DivideAndCapture
