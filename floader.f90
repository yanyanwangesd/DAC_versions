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
  character*400 cread, mstr, mstr2

  ! for each read operation, s holds the return status: 0=OK, else error
  integer ia,s,sf, ci,numpol, nc,n,ib, iok,ic, nset, ncon, readpos, nrast, npolc
  integer hasnan, avgmc, nsetp, nsetr, rformat
  real ra, r1, r2, r3,r4, maxm,minm, avgm
  double precision da, db
  character*400 ca
  character(100000)  buffer
  character(256) filename
  logical filexist
  integer size, stat, pos, rp1, rp2,rpn, pf1, pf2, ncol, nrow, ncells
  integer i, j,k,  filemode
  double precision res, nanrepval,  minx,maxx,miny,maxy,x,y, resx, resy, rrat
  double precision, dimension(:,:), pointer :: matx, xyz




  integer dummy_int
  double precision dummy_double


s=0			! tracks individual read errors (i.e. wrong data type)
mstr=''


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
params%num_bins=iread(s)		                	! how many bins to store erosion values
params%sample_per_bin=iread(s)		              	! how many samples to store in bin

params%read_restart = lread(s)                 	!flag to read GeoFile to restart run from the end of previosu run
params%num_restart = iread(s)			              ! where to restart
if(.not.params%read_restart) params%num_restart = 0 !otherwise it starts the experiment with a read-in step value

params%ascii = lread(s)				                   ! flag to output ASCII data of the run (used for post-processing)
params%show_vtkfine = lread(s)                  ! flag:  output the full triangulation involving the divides. Large output files.

! params%f_varies=.FALSE.
! return

! counted 43 reads so far
! now special provisions
!params%f_varies = lread(s)

! test it
! rewind(10)
! do ia=1,43
!   mstr=cread(s)
! enddo
! params%f_varies = lread(s)
! print*, params%f_varies
! readpos=44;
!
!

  !mstr=cread(s)   ! read f_num_sets

  ! first read pass through input file to determine order/type/existence of input conditions
  ! -------------------------------------------------------------------------
  readpos=43 ! counted 43 reads so far
  s=0
  ic=0
  nrast=0
  npolc=0
  mstr=''
  do while(s.eq.0)
    ci=0
    mstr=cread(s)
    !print*,'s_',s,mstr(1:1),trim(mstr)
    !print*, len_trim(mstr)
    !if ( mstr(1:1).eq.'0' )then
    if ( mstr(1:1).ne.'R' .and. mstr(1:1).ne.'p' )then
       !print*, len_trim(mstr)
       exit
    endif
    readpos=readpos+1
    if ( mstr(1:1).eq.'R' ) then
      ci=1
      do ia=1,6
        mstr=cread(s)
        readpos=readpos+1
      enddo
    else if ( mstr(1:1).eq.'p' ) then
      do ia=1,4
        mstr=cread(s)
        readpos=readpos+1
      enddo
      ia=iread(s)
      do ib=1,ia
        mstr=cread(s)
      enddo
    endif
    if(s.eq.0)then
      ic=ic+1
      if(ci.eq.1)then
        nrast=nrast+1
      else
        npolc=npolc+1
      endif
    endif
  enddo
  print*,'conditions: (tot/poly/raster)',ic,npolc,nrast
  if(ic.gt.0)params%f_varies=.TRUE.
  !stop

  rewind(10)
  do ia=1,43
    mstr=cread(s)
  enddo
  !print*,mstr
  ! -------------------------------------------------------------------------

    params%f_num_sets = ic
    params%f_num_rast = nrast
    params%f_num_polc = npolc

    if(params%f_varies)then

    ! limits: number of conditions specified - none; segments per condition block - 1000; polynomial degree - 50
    allocate(params%f_ctype(params%f_num_sets), params%f_depends_on(params%f_num_sets),params%f_variable_determined(params%f_num_sets), params%f_superpose(params%f_num_sets) )
    allocate(params%f_polyseg(params%f_num_sets))
    allocate(params%pn(params%f_num_sets,1000), params%f_timebound(params%f_num_sets,2), params%f_cidx(params%f_num_sets,2))
    allocate(params%poly(params%f_num_sets,1000,50))
    allocate(params%f_rmarg(params%f_num_rast,4), params%f_rnnode(params%f_num_rast,2))
    allocate(params%f_raster(params%f_num_rast))
    !params%f_raster(params%f_num_rast,4)
    params%f_ctype=0
    params%f_depends_on=0
    params%f_variable_determined=0
    params%f_superpose=0
    params%f_polyseg=0
    params%pn=0
    params%f_timebound=0
    params%f_cidx=0
    params%poly=0
    params%f_rmarg=0
    params%f_rnnode=0
    !params%f_raster=0


    nsetp=0
    nsetr=0
    do nset=1, params%f_num_sets

            !f_input_type
            mstr=cread(s)
            !print*,mstr(1:1)

            !!!!!! this is a polynomial condition !!!!!!!!!!!!!!!!!!!!!!!!!!
            if ( mstr(1:1).eq.'p' ) then
              nsetp=nsetp+1
              params%f_ctype(nset)=1
              print*, '------------  polynomial block',nset,'      -------------'

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
                 print*,trim(mstr)

              enddo
              !!!!!! this is a polynomial condition !!!!!!!!!!!!!!!!!!!!!!!!!!




              !!!!!! this is a raster condition !!!!!!!!!!!!!!!!!!!!!!!!!!
         else
              print*, '------------  raster block',nset,'      -----------------'

              nsetr=nsetr+1
              params%f_ctype(nset)=2
              params%f_cidx(nset,2)=nsetr   ! maps raster index into global condition index

              !filename of raster
              mstr=cread(s)
              !print*,trim(mstr)
              filename=mstr

              filexist=.FALSE.
              inquire(file=filename, exist=filexist)
              !print*, filexist
              if(.not.filexist)then
                  print*, 'ERROR - file does not exist: ',filename
                  stop
              endif

              !determine raster format and read
              rformat=iread(s)
              select case(rformat)
                  case (0)
                    !! --- 0 flat txt matrix
                    open(unit = 11 , file = filename, access='sequential',form='formatted')
                    i=1
                    do while (1.eq.1)
                        read(11,"(A)", ADVANCE='NO', IOSTAT=stat,  SIZE=size, end=999) buffer

                        if(i.eq.1)then
                            rpn=len(trim(buffer))
                            j=0
                            pf1=1
                            rp2=1
                            pos=1
                            rp1=1
                            do while (rp1.lt.rpn)
                                if( buffer(rp1:rp1).eq.' '  .or. buffer(rp1:rp1).eq.',' )then
                                    !skip this one for reading
                                    rp1=rp1+1
                                endif
                                if(buffer(rp1:rp1).ne." "  .and. buffer(rp1:rp1).ne.",")then
                                    !found the start of a value, first non-delim c
                                    if(rp1.eq.1)then
                                        rp2=rp1
                                        j=j+1
                                    else
                                      if(buffer(rp1-1:rp1-1).eq." "  .or. buffer(rp1-1:rp1-1).eq."," )then
                                          rp2=rp1
                                          j=j+1
                                      endif
                                    endif
                                    rp1=rp1+1
                                endif
                            enddo
                        endif

                        i=i+1
                    enddo

                    999 nrow=i-1
                        ncol=j
                        print*, nrow,ncol

                    allocate(params%f_raster(nsetr)%matx(nrow,ncol))
                    params%f_raster(nsetr)%matx=0
                    rewind(11)

                    READ(11,*) ((params%f_raster(nsetr)%matx(i,j),j=1,ncol),i=nrow,1,-1)
                    close(11)
                  !   do i=1,nrow
                  ! !      print*,matx(i,:)
                  !       print*,params%f_raster(1)%p(i,:)
                  !   enddo
                    !! --- 0 flat txt matrix


                    !stop
                  case (1)
                    !! --- 1 asc file
                    open(unit = 11 , file = filename, access='sequential',form='formatted')
                        read(11,"(A)")buffer
                        read(buffer(6:len_trim(buffer)),*)ncol
                        read(11,"(A)")buffer
                        read(buffer(6:len_trim(buffer)),*)nrow
                        read(11,"(A)")buffer
                        read(11,"(A)")buffer
                        read(11,"(A)")buffer
                        read(buffer(9:len_trim(buffer)),*)res
                        read(11,"(A)")buffer

                        !not all asc have nanrepval ...
                        hasnan =0
                        if(buffer(1:6).eq."NODATA" .or. buffer(1:6).eq."nodata")then
                            read(buffer(13:len_trim(buffer)),*)nanrepval
                            print*, nrow, ncol, res, nanrepval
                            hasnan=1
                        else
                            rewind(11)
                             read(11,"(A)")buffer
                            read(11,"(A)")buffer
                            read(11,"(A)")buffer
                            read(11,"(A)")buffer
                            read(11,"(A)")buffer
                            read(buffer(9:len_trim(buffer)),*)res
                            print*, nrow, ncol, res
                        endif

                        allocate(params%f_raster(nsetr)%matx(nrow,ncol))
                        params%f_raster(nsetr)%matx=0
                        READ(11,*) ((params%f_raster(nsetr)%matx(i,j),j=1,ncol),i=nrow,1,-1)
                        close(11)

                        ! now undo the nan values - two loops
                        minm = 1e12
                        maxm = -1e12
                        avgm = 0
                        avgmc = 0
                        do j=1,ncol
                            do i=1,nrow
                                if(params%f_raster(nsetr)%matx(i,j).ne.nanrepval)then
                                    if (  params%f_raster(nsetr)%matx(i,j).gt.maxm )maxm = params%f_raster(nsetr)%matx(i,j)
                                    if (  params%f_raster(nsetr)%matx(i,j).lt.minm )minm = params%f_raster(nsetr)%matx(i,j)
                                    avgm = avgm + params%f_raster(nsetr)%matx(i,j)
                                    avgmc = avgmc + 1
                                endif
                            enddo
                        enddo
                        avgm = avgm / avgmc
                        do j=1,ncol
                            do i=1,nrow
                                if(params%f_raster(nsetr)%matx(i,j).eq.nanrepval)then
                                    params%f_raster(nsetr)%matx(i,j) =avgm
                                endif
                            enddo
                        enddo

                    !! --- 1 asc file

                  case (2)
                    !! --- 2 csv file
                    open(unit = 11 , file = filename, access='sequential',form='formatted')
                        i=0
                        do while (1.eq.1)
                            read(11,"(A)", ADVANCE='NO', IOSTAT=stat,  SIZE=size, end=998) buffer
                            i=i+1
                        enddo

                 998    ncells=i

                        allocate(xyz(ncells,3))
                        rewind(11)
                        READ(11,*) ((xyz(i,j),j=1,3),i=1,ncells)

                        minx=1.d12
                        miny=1.d12
                        maxx=-1.d12
                        maxy=-1.d12
                        do i=1,ncells
                            x=xyz(i,1)
                            y=xyz(i,2)
                            if(x.gt.maxx)maxx=x
                            if(x.lt.minx)minx=x
                            if(y.gt.maxy)maxy=y
                            if(y.lt.miny)miny=y
                        enddo
                        resx=(maxx-minx)
                        resy=(maxy-miny)
                        res=sqrt((resx*resy)/(ncells))
                        ncol=resx/res
                        nrow=ncells/ncol
                        print*, "guessed raster size:",nrow,ncol
                        resx=(maxx-minx)/(ncol-1)
                        resy=(maxy-miny)/(nrow-1)
                        print*,"resolution", resx, resy

                        allocate(params%f_raster(nsetr)%matx(nrow,ncol))
                        !READ(xyz(:,3),*) ((matx(i,j),j=1,ncol),i=1,nrow)
                        k=1
                        do i=1,nrow
                            do j=1,ncol
                                !matx(i,j)=xyz(k,3)
                                params%f_raster(nsetr)%matx(i,j)=xyz(k,3)
                                k=k+1
                            enddo
                        enddo
                        deallocate(xyz)
                        close(11)
                    !! --- csv file
!                    print*, 'csv reader not yet implemented'
!                    stop
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

              params%f_rnnode(nsetr,1)=nrow
              params%f_rnnode(nsetr,2)=ncol


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

              ! margin offlaps
              mstr=cread(s)
              read (mstr,*,iostat=iok) params%f_rmarg(nsetr,:)!r1,r2,r3,r4
              print*, 'margins', params%f_rmarg(nsetr,:)!r1,r2,r3,r4


              !stop
              !!!!!! this is a raster condition !!!!!!!!!!!!!!!!!!!!!!!!!!
         endif ! if polynome/raster
    enddo ! loop all blocks
    ! f_num_set - blocks
    ! polyseg - lines
    ! pn    -  entries in line

     print*, '----------------- ',params%f_num_sets,'     blocks set  ----------------'

endif !(if_params%f_varies)

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

	character*400 line
	double precision da
	integer return_status, rs, ffread

	return_status=ffread(line)

	read(unit=line,fmt=913,err=713) da
913	format(D400.140)

	dread=da

	return

713	print*,'FAILED to read double, got:   ',line
	rs=1
  stop

end function dread




!----------------------------------------------------------
real function rread(rs)

implicit none

	character*400 line
	real ra
	integer return_status, rs, ffread

	return_status=ffread(line)

	read(unit=line,fmt=914,err=714) ra
914	format(E400.140)

	rread=ra

	return

714	print*,'FAILED to read real, got:   ',line
	rs=1
  stop

end function rread




!----------------------------------------------------------
integer function iread(rs)

implicit none

	character*400 line
	integer ia
	integer return_status, rs, ffread

	return_status=ffread(line)

	!handle possible error
	if(return_status.ne.0) then
		rs=1
		return
	end if


	read(unit=line,fmt=915,err=715) ia
915	format(I400.140)

	iread=ia

	return

715	print*, 'FAILED to read integer, got:   ',line
	rs=1
  stop

end function iread




!----------------------------------------------------------
character*400 function cread(rs)

implicit none

	character*400 line
	character*400 cline
	integer return_status, rs, ffread, iok

  iok=0

	return_status=ffread(line)

  if(return_status.ne.0) return

	read(unit=line,fmt=916,err=716,iostat=iok) cline
916	format(A400)


 if (iok.ne.0)then
   rs=1
 else
	cread=cline
  endif
	return

716	print*, 'FAILED to read string, got:   ',line
	rs=1
  stop


end function cread





!----------------------------------------------------------
logical function lread(rs)

implicit none

	character*400 line
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
  stop

end function lread





!----------------------------------------------------------
! master function, actual read from file

integer function ffread(line)

implicit none

	character*400 line
	character*400 fline
	character firstc
	integer i
	integer state, ios

	state=0
  ios=0

	do while(state.eq.0)
    !read (unit=10, fmt=911,err=711, iostat=ios, end=789) fline
    read (unit=10, fmt=911, iostat=ios, end=789) fline
911		format(A400)
    if (ios.lt.0) then
      print*,'offensive content:', fline
      !return
    endif

		read(fline,912, err=709, end=789)firstc
912		format(A1)
913   format(A)

		!if (firstc.eq.'/' .or. firstc.eq.' ' .or. firstc.eq.achar(10) .or. firstc.eq.new_line(firstc)) then
    if (firstc.eq.'/' .or. firstc.eq.' ' ) then
			!print*, 'encountered a comment line'
		else
			state=1
			line = fline
		endif
	end do

	ffread=0

	!print*, 'done ff'
     !print*, line

	return

709	ffread=1
      print*, 'failed low-level read from string - exiting'
      stop

711	ffread=1
    print*, 'low-level read operation failed - exiting'
    print*, trim(fline)
    stop
	! return an error

789 ffread=2!print*, ' '!'reached end of file'


end function ffread

! End of modified section ------------------
