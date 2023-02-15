subroutine floader(  params, geometry )

use definitions

implicit none

 type (parm) params
 type (geom) geometry
  type (del) delaunay
  type (stck) stack
  type (netw) network

  !function "prototypes"
  integer iread
  double precision dread
  real rread
  logical lread
  character*150 cread, mstr, mstr2

  ! for each read operation, s holds the return status: 0=OK, else error
  integer ia,s,sf, ci,numpol, nc,n,ib, iok,ic, nset
  real ra, r1, r2, r3,r4
  double precision da, db
  character*150 ca


  
  integer dummy_int
  double precision dummy_double
  

s=0			! tracks individual read errors (i.e. wrong data type)
  
  
  ! CRUCIAL  -  success of opening decides which section of code is executed	
  open (unit=10, file='init.dac', status='old', iostat=sf, err=701)
  
  
    ! do all the loading, sequentially stored in init.dac



 
geometry%xl=dread(s)
geometry%yl=dread(s)
geometry%nx=iread(s)
geometry%ny=iread(s)
geometry%zl=dread(s)
geometry%nnmax=iread(s)
params%xc=dread(s)                             	! hillslope length (in the vicinity of divides) (in m)
params%lmax=dread(s)                           	! maximum distance between two points before adding a node
params%amin=dread(s)                           	! minimum catchment head area before node remove_d
params%ldivmax=dread(s)                        	! maximum divide length for adding nodes
params%max_adv = dread(s)                      	! max advection velocity
do ia=1,4
  geometry%boundary(ia)=iread(s)
enddo

 	  
params%deltat=dread(s)                        	! time step length (in years)
params%time=dread(s)                           	! time at start (should always be 0)                  
params%tfinal=dread(s)                        	! final time (in years)
params%freq=iread(s)                     	! frequency of plots in number of time steps 
params%istep=0                                  ! time step counter (should always be set to 0)
params%k_scalar1=dread(s)                      	! erodibility constant
params%n=dread(s)                              	! slope exponent in fluvial incision law
params%m=dread(s)                             	! area (or discharge) exponent in fluvial incision law
params%h=dread(s)                            	! Hack's exponent
params%hmn=1.d0-params%h*params%m/params%n     	! exponent combination parameter (do not change)
params%ka=dread(s)                            	! Hack's constant (no unit)
params%tanthetac=dread(s)                   	! tangent of hillslope slope
params%diffusivity=dread(s)                    	! diffusivity for hilltops
params%diffusivity=params%diffusivity*2.d0     	! double diffusivity as it always appears with a 2*
params%min_erosion_rate=dread(s)               	! minimum allowable erosion rate for diffusion channel head
params%min_tan_head_slope = dread(s)         	! minimum slope of channel head (used as a threshold to diffusion) 
params%rainfall_height=dread(s)              	! this is used for constant rainfall parameter
params%uplift_scalar1=dread(s)                 	! uplift rate as constant in interior
params%uplift_scalar2=params%uplift_scalar1     ! uplift rate on boundary (used only for channel head calc)

params%add_nodes=lread(s)                      	! flag to allow adding nodes
params%move_points=lread(s)                   	! flag to allow the horizontal motion of nodes
params%capture=lread(s)                       	! flag to allow for captures
params%divide=lread(s)                        	! flag to allow for divide calculations
params%small_divide=lread(s)                  	! flag to allow for small divide calculations (case where divide is very close to zj)
params%diffusion=lread(s)                      	! flag to allow diffusive hilltops
params%transient_divide=lread(s)               	! flag to calculate transient elevation of fluvial part of divides
params%num_bins=iread(s)			! how many bins to store erosion values
params%sample_per_bin=iread(s)			! how many samples to store in bin

params%read_restart = lread(s)                 	!flag to read GeoFile to restart run from the end of previosu run
params%num_restart = iread(s)			! where to restart

params%ascii = lread(s)				! flag to output ASCII data of the run (used for post-processing)
params%show_vtkfine = lread(s)                  ! flag:  output the full triangulation involving the divides. Large output files.





! now special provisions
params%f_varies_with_xyz = lread(s)

if(params%f_varies_with_xyz)then

    params%f_num_sets = iread(s)

    ! limits: number of conditions specified - none; segments per condition block - 1000; polynomial degree - 50
    allocate(params%f_polyseg(params%f_num_sets), params%f_depends_on(params%f_num_sets),params%f_variable_determined(params%f_num_sets), params%f_superpose(params%f_num_sets) )
    allocate(params%pn(params%f_num_sets,1000), params%f_timebound(params%f_num_sets,2))
    allocate(params%poly(params%f_num_sets,1000,50))




    do nset=1, params%f_num_sets

            print*, '-----------------  block',nset,'      ----------------'

            ! f_depends_on          0-x, 1-y, 2-z
            mstr=cread(s)
            select case(trim(mstr))
                case ("x")
                    print*, 'independent variable x'
                    params%f_depends_on(nset)=0
                case ("y")
                    print*, 'independent variable y'
                    params%f_depends_on(nset)=1
                case ("z")
                    print*, 'independent variable z'
                    params%f_depends_on(nset)=2
                case default
                    print*, "ERROR - unknown independent variable case: ",trim(mstr), ',  block:', nset
                    stop
            end select
            ! f_variable_determined 0-u, 1-v, 2-w, 3-P, 4-k
            mstr=cread(s)
            select case(trim(mstr))
                case ("u")
                    print*, 'dependent variable u'
                    params%f_variable_determined(nset)=0
                case ("v")
                    print*, 'dependent variable v'
                    params%f_variable_determined(nset)=1
                case ("w")
                    print*, 'dependent variable w'
                    params%f_variable_determined(nset)=2
                case ("p" , "P")
                    print*, 'dependent variable p'
                    params%f_variable_determined(nset)=3
                case ("k")
                    print*, 'dependent variable k'
                    params%f_variable_determined(nset)=4
                case ("z")
                    print*, 'dependent variable z'
                    if(params%f_depends_on(nset).eq.2) stop 'cannot set z as function of itself'
                    params%f_variable_determined(nset)=5
                 case default
                    print*, "ERROR - unknown dependent variable case: ",trim(mstr), ',  block:', nset
                    stop
            end select

            ! how to combine overlapping conditions
            params%f_superpose(nset) = iread(s)

            ! timebounds
            mstr=cread(s)
            ib=len(trim(mstr))
            ic=1
            do while(mstr(ic:ic).ne.',' .and. ic.lt.149)
              ic=ic+1
            enddo
            if(ic-1.eq.ib  .or. ic.eq.149)stop "wrong input for time"
            read (mstr(1:ic-1),*,iostat=iok) da
            read (mstr(ic+1:ib),*,iostat=iok) db
            print*, 'times: ', da,db
            params%f_timebound(nset,1)=da
            params%f_timebound(nset,2)=db

            ! how many segments
            numpol = iread(s)
            if(numpol.gt.1000)stop 'currently only up to 1000 segments supported'
            params%f_polyseg(nset) = numpol  !+1         ! number of piecewise segments

            ! read segment polynomial coordinates
            do ia=1, numpol
               mstr=cread(s)

               ib=len(trim(mstr))
               n=1
               nc=1
               do ic=1,ib
                  if(mstr(ic:ic).eq.',' .or. ic.eq.ib)then
                      read (mstr(nc:ic-1),*,iostat=iok) da
                      if(ic.eq.ib)read (mstr(nc:ic),*,iostat=iok) da
                      params%poly(nset, ia,n)=da
                      n=n+1
                      nc=ic+1
                  endif
               enddo
               params%pn(nset, ia)=n-1
               print*,mstr

            enddo
    enddo
    ! f_num_set - blocks
    ! polyseg - lines
    ! pn    -  entries in line

     print*, '----------------- ',params%f_num_sets,'     blocks set  ----------------'

endif

close(10)

!stop
return

!	error on file-existance level
701	if(sf.ne.0) then
		print*, 'WARNING : No input file found, or reading failed'
		print*,''
		print*,'awaiting acknowledgement - hit return'
		read(*,*)
	end if
	  

end subroutine floader
























!----------------------------------------------------------
double precision function dread(rs)
	
implicit none

	character*150 line
	double precision da
	integer return_status, rs, ffread
	
	return_status=ffread(line)
	
	read(unit=line,fmt=913,err=713) da
913	format(D150.140)
	
	dread=da

	return

713	print*,'FAILED to read double, got:   ',line
	rs=1
	
end function dread




!----------------------------------------------------------
real function rread(rs)
	
implicit none

	character*150 line
	real ra
	integer return_status, rs, ffread
	
	return_status=ffread(line)
	
	read(unit=line,fmt=914,err=714) ra
914	format(E150.140)
	
	rread=ra

	return
	
714	print*,'FAILED to read real, got:   ',line
	rs=1
	
end function rread




!----------------------------------------------------------
integer function iread(rs)
	
implicit none

	character*150 line
	integer ia
	integer return_status, rs, ffread
	
	return_status=ffread(line)
	
	!handle possible error
	if(return_status.ne.0) then
		rs=1
		return
	end if
	
	
	read(unit=line,fmt=915,err=715) ia
915	format(I150.140)
	
	iread=ia

	return
	
715	print*, 'FAILED to read integer, got:   ',line
	rs=1
	
end function iread




!----------------------------------------------------------
character*150 function cread(rs)
	
implicit none

	character*150 line
	character*150 cline
	integer return_status, rs, ffread, iok
	
	return_status=ffread(line)
	
	read(unit=line,fmt=916,err=716,iostat=iok) cline
916	format(A150)
	
	cread=cline

	return
	
716	print*, 'FAILED to read string, got:   ',line
	rs=1
	
end function cread





!----------------------------------------------------------
logical function lread(rs)
	
implicit none

	character*150 line
	integer ia
	logical la
	integer return_status, rs, ffread
	
	return_status=ffread(line)
	
	
	
	
	! VERSION 1 - read logical in FORTRAN language (too awkward for DAC input)
	!read(unit=line,fmt=917,err=717) la
!917	format(L10)	
	!lread=la
	
		
		
	! VERSION 2 - read logical as 0;1
	read(unit=line,fmt=917,err=717) ia
917	format(I30)
	la=.TRUE.
	if(ia.eq.0) la=.FALSE.
	if(ia.eq.-1)la=.TRUE.
	lread=la

	
	
	return
	
717	print*,'FAILED to read logical, got:   ',line
	rs=1
	
end function lread





!----------------------------------------------------------
! master function, actual read from file
	
integer function ffread(line)

implicit none

	character*150 line
	character*150 fline
	character firstc
	integer i
	integer state
		
	state=0
	
	do while(state.eq.0)
		read (unit=10, fmt=911,err=711) fline
911		format(A150)
		read(fline,912)firstc
912		format(A1)
		
		if (firstc.eq.'/' .or. firstc.eq.' ') then
			!print*, 'encountered a comment line'
		else
			state=1
			line = fline	
		endif
	end do
	
	ffread=0
	
	!print*, 'done ff'
	
	return

711	ffread=1
	! return an error


end function ffread
	
! End of modified section ------------------	
