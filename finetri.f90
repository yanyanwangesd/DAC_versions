subroutine finetri (geometry, params,delaunay, network, landscape)


  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay
  type (netw) network
  type (stck) stack
  type (parm) params
  TYPE (lands) landscape
  


!! declarations  
  integer i,j,k,l, m,temp1,temp2
  integer nd, nfp, nlocp,di(3,2),   divno,  trialg, divstyle
  integer tix(3), rix(3), fix(3),divlocix								
  double precision xd(3), yd(3),zd(3), ctrx,ctry,ctrz
  double precision dist(3), disttot

  real*8,dimension(:,:),allocatable::points
  real*8,dimension(:,:),allocatable::centers
  integer,dimension(:,:),allocatable::e,v
  integer,dimension(:),allocatable::vis_tlist,vis_elist,add_tlist
  integer num,numtri_max,numtri,nv_max,mode,nfirst,itstart
  logical,dimension(:),allocatable::inactive,subset
  real*8 eps
  real*8  ctr(2)




  !!  USER CHOICES  

  !! ------------------------------------------------------------------------------
  !!   choose here whether to use delaun or simple own triangulation
	trialg = 0
	!! trialg 0   ->   use simple fix triangle algorithm
	!! trialg 1   ->   use delaun in library
  !! ------------------------------------------------------------------------------


  !! ------------------------------------------------------------------------------
  !!   choose here whether to use divides as calculated or straight in case of two divides
	 divstyle = 0
	!! divstyle 0   ->   use straight connection for two divides
	!! divstyle 1   ->   use centroid
  !! ------------------------------------------------------------------------------







!! set-up for call to delaun in libnn
num=9		! 3 nodes, 1 divide point within, max 3 intersections with divides; twice the chance for boundary triangles
numtri_max=9
nv_max=7
mode=0
nfirst=1
itstart=1
eps=1.d-6



!! local arrays to feed the call to delaun
allocate (points(2,num))
allocate (centers(3,numtri_max))
allocate (e(3,numtri_max),v(3,numtri_max))
allocate (vis_tlist(nv_max),vis_elist(nv_max),add_tlist(nv_max))
allocate (inactive(num),subset(num))


if (landscape%nfinetri.ne.0) deallocate (landscape%xx,landscape%yy,landscape%zz,landscape%edot,landscape%fineicon)

! print*, 'I think that need only', delaunay%ntriangles*10 
allocate (landscape%xx(delaunay%ntriangles*10),landscape%yy(delaunay%ntriangles*10))
allocate (landscape%zz(delaunay%ntriangles*10),landscape%edot(delaunay%ntriangles*10))
allocate(landscape%fineicon(3,delaunay%ntriangles*10))

 landscape%edot=0.
!! calculations --------------------------------------------------------------------------------------------------------------------------------------------------------


landscape%nfinenode = 0
landscape%nfinetri = 0
landscape%xx=0
landscape%yy =0
landscape%zz=0

!! loop triangles
do i = 1, delaunay%ntriangles
	! how many divides per given triangle: nd
	nd = delaunay%numdivides(1,i)
	num = 3 + 1 + nd
	nfp = num
	

	nlocp = 0
	points = 0
	di = 0
	xd = 0
	yd = 0
	zd = 0
	divno = 0
	ctr = 0
	tix = 0
	rix = 0
	
	
	! load triangle vertices into local point array
	tix(:) = delaunay%icon(1:3,i) 
	points(1,1:3) = geometry%x(  tix )
	points(2,1:3) = geometry%y(  tix  )
   

	
	
	! get divide end point indices (di) and coordinates (xd,yd,zd)
	nlocp = 4    !! three edges stored, one center to come
	do j=1, nd		! loop divides
		
		divno = delaunay%numdivides(  j+1,i  )
		di(j,1) = geometry%nndivnode(  divno ,1   )
		di(j,2) = geometry%nndivnode(  divno, 2   )

		xd(j) = geometry%xdiv(  divno  )
		yd(j) = geometry%ydiv(  divno  )
		zd(j) = geometry%zdiv(  divno  )
		
		points(1,nlocp+1) = xd(j)
		points(2,nlocp+1) = yd(j)
		
		nlocp = nlocp+1
		
	enddo
	
	
	!! Chosen mid-point
	!! add center point to point array, update fill state of array,
	if(  divstyle.eq.1 ) then
		call calculate_centroid (points, ctr)
	endif
	
	!!  STRAIGHT 2-divide connections -  center = mean of divide points
	if(  divstyle.eq.0) then
		if (nd.eq.2 )then
			ctr(1) = ( points(1 ,5) + points(1,6) ) /2
			ctr(2) = ( points(2 ,5) + points(2,6) ) /2
			
		else
			call calculate_centroid (points, ctr)
		endif
	endif
	
	
	
	
	points(1,4) = ctr(1)
	points(2,4) = ctr(2)

	

	
	num = nlocp
	nfp = nlocp
	
	
	!! triangulation is done in 2d only, so need to obtain z-values
	ctrz=0
	disttot = 0
	
	if(nd.eq.1) THEN
		ctrz = zd(1)
	ELSE
		ctrx = ctr(1)
		ctry = ctr(2)
		!! weighted mean
		DO j = 1,nd
			dist(j) = dsqrt(      (ctrx-xd(j))**2      +       (ctry-yd(j))**2       )
			disttot = disttot + 1/dist(j)
			ctrz = ctrz+zd(j)/dist(j)
		ENDDO
			ctrz = ctrz/disttot
	ENDIF

	!! in rare case there is no divide
	if(nd.eq.0) ctrz = (geometry%z(delaunay%icon(1,i))+geometry%z(delaunay%icon(2,i))+geometry%z(delaunay%icon(3,i)))/3

	!! sanitize
	if(ctrz.gt.10000) ctrz =0
	


	!! local triangulation from points

	if (trialg .eq.1) then
		call delaun (points,num,e,v,numtri,numtri_max, &
		     vis_tlist,vis_elist,add_tlist,eps,nv_max, &
		     mode,inactive,nfirst,itstart,subset)
	endif
	
	if(trialg .eq.0) then
		
		! use simple custom triangulation
		numtri = 0
		!first check whether nodes are river nodes
		do j=1,3
			k = mod(j,3)+1        !! 2,3,1 
			! check for each node index whether it has a river connection with its clockwise adjacent corner --  GLOBAL INDICES tix,network%receiver, LOCAL index rix for sides 
			if(     ( network%receiver( tix(k) ) .eq. tix(j) ) .or. (network%receiver( tix(j) ) .eq. tix(k) )   ) then
				rix(j) = 1
			endif
			! store which ones are boundaries
			if(geometry%fix(tix(j)) .eq. 1) fix(j) = 1
		end do
		
		! calculate triangles depending on configuration
		select case (nd)
			case(0)		!! connect all triangle nodes with center  (case  0 happens on weird occasions )
				do j=1,3
					k = mod(j,3)+1    !! 2,3,1      4 = centroid
					v(1:3, j) = (/  j,k, 4 /) 
					numtri = 3
				enddo
				
				
			case(1)		!! connect triangles with centers if they are river nodes 
				do j=1,3
					numtri = numtri +1
					k = mod(j,3)+1    !! 2,3,1      4 = centroid, 5 = single divide
					if (rix(j).eq.1)then	! river-bound side, between node j and j+1 (circular indexing)
						v(1:3, numtri) = (/  j,k, 4 /)
					else
						v(1:3, numtri) = (/ j,5, 4  /)
						numtri = numtri +1
						v(1:3, numtri) = (/ 5,k,4  /)
					endif
				enddo
				
				!! treat boundary, connect two bdy nodes to ctr
				If ( sum(fix).gt.0 ) then
						if(     fix(1).eq.1 .and. fix(2).eq.1     )   then
							numtri = numtri +1
							v(1:3, numtri) = (/ 1,2,4  /)
						end if
						If (     fix(1).eq.1 .and. fix(3).eq.1     )  THEN
							numtri = numtri +1
							v(1:3, numtri) = (/ 1,3,4  /)
						end if
						if (     fix(2).eq.1 .and. fix(3).eq.1     ) THEN
							numtri = numtri +1
							v(1:3, numtri) = (/ 2,3,4  /)
						endif
				endif
				
				
				
			case(2)		!! connect river side with center, nodes with clockwise a non-river next  with divide, center, and neighbour
				do j=1,3
					numtri = numtri +1
					k = mod(j,3)+1    !! 2,3,1      4 = centroid, 5,6 divides
					if (rix(j).eq.1)then	! river-bound side, between node j and j+1 (circular indexing)
						v(1:3, numtri) = (/  j, k,4 /)
					else
						! which one of both divides -look in nndivnodes
						l = delaunay%numdivides(  2,i  )	! if it is not this one, then the other = delaunay%numdivides(  3,i  )
						divlocix = 1
						temp1 = geometry%nndivnode(l,1)
						temp2 = geometry%nndivnode(l,2)
if(.not.(((temp1.eq.tix(j)).and.(temp2.eq.tix(k))).or.((temp1.eq.tix(k)).and.(temp2.eq.tix(j)))))then
							l = delaunay%numdivides(  3,i  )
							divlocix = 2
						endif
						
						v(1:3, numtri) = (/ j,4+divlocix, 4  /)
						numtri = numtri +1
						v(1:3, numtri) = (/ 4+divlocix,k,4  /)
					endif
				enddo
			
				!! treat boundary, connect two bdy nodes to ctr
				If ( sum(fix).gt.0 ) then
						if(     fix(1).eq.1 .and. fix(2).eq.1     )   then
							numtri = numtri +1
							v(1:3, numtri) = (/ 1,2,4  /)
						endif
						if (     fix(1).eq.1 .and. fix(3).eq.1     )  THEN
							numtri = numtri +1
							v(1:3, numtri) = (/ 1,3,4  /)
						 endif
						 if (     fix(2).eq.1 .and. fix(3).eq.1     ) THEN
							numtri = numtri +1
							v(1:3, numtri) = (/ 2,3,4  /)
						endif
				endif
				
				
			case(3)		!! connect all nodes and divides with centers 
				do j=1,3
					numtri = numtri +1
					k = mod(j,3)+1    !! 2,3,1      4 = centroid, 5,6,7 divides
					! which one of the divides is between the node and its clockwise neighbour -look in nndivnodes
					do m=1,3 
					l = delaunay%numdivides(m+1,i)
					temp1 =geometry%nndivnode(l,1)
					temp2 =geometry%nndivnode(l,2)
					if (((temp1.eq.tix(j)).and.(temp2.eq.tix(k))).or.((temp1.eq.tix(k)).and.(temp2.eq.tix(j)))) then
						divlocix = m
					endif
					enddo
					
					v(1:3, numtri) = (/ j,4+divlocix, 4  /)
					numtri = numtri +1
					v(1:3, numtri) = (/ 4+divlocix,k,4  /)
				enddo
				
				
		end select
		
	endif




	!! update from delaun routine
	! update triangle-to-node connectivity, shift to match chunk   = fill new  "icon" array
	landscape%fineicon(1:3, landscape%nfinetri + 1 : landscape%nfinetri + numtri ) = v(:,1:numtri) + landscape%nfinenode 
	landscape%nfinetri = landscape%nfinetri + numtri

	! fill new "node" array
	landscape%xx (  landscape%nfinenode+1 : landscape%nfinenode+nlocp ) = points(1, 1:nlocp)
	landscape%yy (  landscape%nfinenode+1 : landscape%nfinenode+nlocp ) = points(2, 1:nlocp)
	! z not explicitly used before, so piecewise
	landscape%zz (  landscape%nfinenode+1 : landscape%nfinenode+3 ) = geometry%z(  delaunay%icon(1:3,i)  )
	landscape%zz ( landscape%nfinenode+4  ) = ctrz
	landscape%zz (  landscape%nfinenode+5 : landscape%nfinenode+nlocp ) = zd(1:nd)
       ! erosion rate not explicitly used before, so piecewise
        landscape%edot (  landscape%nfinenode+1 : landscape%nfinenode+3 ) = geometry%erosion_rate(  delaunay%icon(1:3,i)  )
        !print*, geometry%erosion_rate(  delaunay%icon(1:3,i)  ), 'tri'
        landscape%edot ( landscape%nfinenode+4  ) = SUM(geometry%erosion_rate(  delaunay%icon(1:3,i)  ))/3.
        !print*, SUM(geometry%erosion_rate(  delaunay%icon(1:3,i)  ))/3., 'center'
        do j=1,nd
           landscape%edot(landscape%nfinenode+j+4)=(geometry%erosion_rate(di(j,1))+geometry%erosion_rate(di(j,2)))/2.
           !print*, (geometry%erosion_rate(di(j,1))+geometry%erosion_rate(di(j,2)))/2.,'div'
        enddo
	landscape%nfinenode = landscape%nfinenode + nlocp
	


enddo

!print*, 'Delaunay tri=', delaunay%ntriangles, 'Fine tri=',landscape%nfinetri


!! deallocate local arrays
deallocate (points,centers,e,v,vis_tlist,vis_elist,add_tlist,inactive,subset)



return

end subroutine finetri

