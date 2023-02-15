MODULE definitions

   type parm
      integer istep,freq,num_restart
      double precision deltat,time,tfinal,h,ka,xc,tanthetac,hmn,m,n
      double precision rainfall_available,rainfall_minimum,rainfall_height
      double precision lmax,diffusivity,amin,min_erosion_rate,uplift_scalar1
      double precision uplift_scalar2,k_scalar1,ldivmax, max_adv, min_tan_head_slope
      logical plot_triangles,plot_receiver,plot_donors,plot_no_receiver
      logical plot_no_donors,write_discharge,write_stack_order,plot_catchment
      logical plot_height,plot_precipitation
      logical move_points,capture,divide,small_divide,diffusion
      logical add_nodes
      logical read_restart, ascii, transient_divide
      logical show_vtkfine
      integer num_bins, sample_per_bin
      logical f_varies_with_xyz
      integer f_num_sets        ! number of condition blocks
      integer,dimension(:),pointer:: f_depends_on, f_variable_determined,f_polyseg, f_superpose     ! dimension f_num_sets
      double precision,dimension(:,:),pointer:: f_timebound       ! limit validity of condition block in time
      double precision,dimension(:,:,:),pointer::poly       ! polynomial coefficients
      integer,dimension(:,:),pointer:: pn   ! degree/elements per line
   end type parm

   type geom
      double precision xl,yl,zl
      double precision,dimension(:),pointer::x,y,z,u,v,w
      double precision,dimension(:),pointer::xdiv,ydiv,zdiv
      double precision,dimension(:),pointer::surface,discharge,precipitation
      double precision,dimension(:),pointer::k,erosion_rate,sediment_flux, chi
      integer,dimension(:),pointer::strahler
      double precision,dimension(:,:),pointer::surface_share
      integer,dimension(:),pointer::fix,nb,catchment
      integer,dimension(:,:),pointer::nn
      integer,dimension(:,:),pointer::nndivtri !for each divide two triangles  
      integer,dimension(:,:),pointer::nndivnode !for each divide two nodes
      integer,dimension(:),pointer::nt ! number of triangles for each node
      integer, dimension(:,:),pointer::nnodetri ! for each node its neighboring triangles
      double precision, dimension(:,:),pointer::erosion_rate_history
      integer nx,ny,nnode,nnode_max,nnmax,ndivide,ndivide_max,ncapture
      integer boundary(4)
   end type geom

   type del
      integer ntriangles
      integer,dimension(:,:),pointer::icon,neighbours,numdivides
      double precision,dimension(:,:),pointer::centers
   end type del

   type netw
      integer nnode,ndonmax,nlake
      integer,dimension(:),pointer::receiver,ndon,lakes,lakes_catch
      integer,dimension(:,:),pointer::donors
   end type netw

   type stck
      integer nnode
      integer,dimension(:),allocatable::order
   end type stck

   type timr
      sequence
      integer ntimer
      character*256 name(1024)
      real time_spent(1024)
      character*256 namein
      real timein
      real time_start,time_stop
   end type timr

   TYPE lands
	sequence
	integer nfinetri, nfinenode
	double precision, dimension(:),allocatable::xx,yy,zz,edot
	integer, dimension(:,:), pointer::fineicon
   END TYPE lands

end MODULE definitions
