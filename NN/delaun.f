c
c------------------------------------------------------------------------
c
c	delaun - calculates delaunay triangulation incrementally 
c	 	 for a set of points in 2-D using a variation of
c		 Lawson's algorithm.
c
c	Input:
c		points(2,np)		array of node co-ordinates
c		num			number of nodes to be used
c               vis_tlist(nv_max)       List of triangles visible from 
c					current point.
c               vis_elist(nv_max)       List of edges visible from 
c					current point.
c               add_tlist(nv_max)       work array used by routine addpoint
c               eps                     distance from an interface for a
c                                       a point to be considered on an
c                                       interface (real*8). Prevents zero
c                                       area triangles resulting from rounding
c                                       error when nodes are co-linear.
c		nv_max			size of work arrays
c		mode			(=0,1,2,3) operation mode (see below)
c		inactive(np)		logical array. If mode=1 then the i-th
c					node is ignored if active(i) = .true.
c		nfirst			If mode=3 then nfirst is the first
c					node added to an existing triangulation
c		itstart			If mode=3 then itstart is a first
c					guess triangle containing node first node
c		subset(np)		logical array. If mode=2 then only
c					the nodes (subset(i),i=1,num) are used.
c
c	Output:
c               v(3,*)           	array of triangle vertices
c               numtri                  number of triangles in current
c                                       triangulation.
c               e(3,*)                  adjacency matrix of neighbouring
c                                       triangles. e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2) 
c					in triangle j, stored counterclockwise 
c					about j.  
c                                       (This is the `opposite' definition)
c
c	Comments:
c
c       This routine calculates the Delaunay triangulation of a set of nodes 
c	using a variation of Lawson's method.  Each node is added sequentially 
c	and the Delaunay triangulation is updated. If the new node is inside 
c	the convex hull of the existing triangulation then the standard Lawson 
c	method is used. If it is outside then the list of triangle edges 
c	which are visible from the new point is calculated using routine 
c	visiblelist and each of these is used as the start of the swapping 
c	routine addpoint.
c
c	Four different operation modes are allowed.
c
c	MODE = 0:
c	The `standard' mode. All nodes from 1 to num are included. The arrays 
c	`subset' and `inactive' are unused and may be set to dummy variables
c       (saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 1:
c	All nodes from 1 to num are included except those for which
c	inactive(i) is set to true. The array `subset' is unused and may
c	be set to a dummy variable (saving memory). The variables nfirst 
c	and itstart are also unused.
c
c	MODE = 2:
c	Only nodes from subset(1) to subset(num) are included. The array
c	`inactive' is unused and may be set to a dummy variable
c	(saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 3:
c	Used to add nodes from nfirst to num to an existing triangulation.
c	Nodes for which inactive(i) is set to true are ignored.
c	The array `subset' is unused and may be set to a dummy variable.
c
c       The performance may be sensitive to the order in which the nodes are
c       added so these can be sorted before calling this routine if desired.
c
c	This routine was converted to use the `opposite' definition of
c	the adjacency matrix on 30/1/96.
c
c
c	Calls are made to triloc_del,visiblelist,insert_point,addpoint.
c
c					         M. Sambridge, Dec. 1994.
c					Modified by J. Braun, Sept. 1995.
c					(last change 30/1/96: multiple modes,
c					 and uses opposite definition of
c					 adjacency matrix)
c
c------------------------------------------------------------------------
c
	Subroutine delaun (points,num,e,v,numtri,numtri_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     mode,inactive,nfirst,itstart,subset)

	real*8		points(2,*)
	real*8		x,y
        real*8          eps,del1,del2,del
 	integer		vis_tlist(*),vis_elist(*),add_tlist(*)
 	integer		v(3,*)
 	integer		e(3,*)
 	integer		subset(*)
	integer		t,p,ccw
        logical         out
	logical		newpoint
	logical		inactive(*)

        if (mode.eq.0.or.mode.eq.1.or.mode.eq.2) then

c					We are calculating Delaunay 
c					of all input points or a
c					subset of all input points

c					find first two active nodes
           if(mode.eq.0)then
              i1=1 
              i2=2 
              nodestart=3
           else if(mode.eq.1)then
              i1=0
              i2=0
              do i=1,num
                 if (i2.ne.0) goto 2222
                 if (i1.ne.0.and..not.inactive(i)) i2=i
                 if (i1.eq.0.and..not.inactive(i)) i1=i
              end do
 2222         continue
              nodestart = i2+1
           else if(mode.eq.2)then
              i1 = subset(1)
              i2 = subset(2)
              nodestart = 3
           end if
c                                       Find three non-colinear points
c                                       to form the first triangle 
           v(1,1) = i1
           v(2,1) = i2
           do 10 j=nodestart,num
              i = j
              if(mode.eq.2)then
                  i = subset(j)
              else if(mode.eq.1.and.inactive(i))then
                  go to 10
              end if
              istart=i
	      del1 = (points(2,i1)-points(2,i))
     &              *(points(1,i2)-points(1,i))
	      del2 = (points(1,i1)-points(1,i))
     &              *(points(2,i2)-points(2,i))
              del = del1-del2
              if(dabs(del).gt.eps) goto 11111
 10        continue
           stop 'all input data are in a line...'
11111      v(3,1) = istart

c					Initialize adjacency matrix
 	   e(1,1) = 0
 	   e(2,1) = 0
 	   e(3,1) = 0
c					Ensure initial triangle 
c					is in ccw order
c					
 	   if(ccw(points(1,v(1,1)),
     &            points(1,v(2,1)),
     &            points(1,v(3,1)),k).eq.-1)then
                  itemp = v(1,1)
                  v(1,1) = v(2,1)
                  v(2,1) = itemp
c	          write(*,*)' initial triangle was cw'
 	    end if
	
c					Initialize variables
 	    numtri = 1
	    t = 1

        else if (mode.eq.3) then
c					We are adding nodes to an 
c					existing triangulation
c					Perform initialization
           nodestart=nfirst
           t = itstart
           istart = 0
           if(t.le.0.or.t.gt.numtri)t=1

        end if 
c					Incrementally update the 
c					Delaunay triangulation

 	do 100 j=nodestart,num


           p = j
           if(mode.eq.1.or.mode.eq.3)then
             if (inactive(j)) goto 100
           else if(mode.eq.2)then
	     p = subset(j)
           end if
             
           if(p.eq.istart)go to 100

	   x = points(1,p)
	   y = points(2,p)

c					locate triangle 
c					containing current node

	   call triloc_del(x,y,points,v,e,t,eps,out,ipos,iface)


	   if(out)then

c					point is outside of convex hull, 
c					so find list of edges that are 
c					visible from current point 

	      call visiblelist(points,e,v,x,y,t,ipos,eps,
     &                         vis_tlist,vis_elist,nvis)

c					for each visible edge 
c					start swapping algorithm

              newpoint = .true.

	      if(nvis.gt.nv_max)then
                 write(*,*)' Error in subroutine delaun:'
                 write(*,*)' Too many visible triangles
     &                       from current point'
                 write(*,*)' Remedy: increase size of parameter nv_max'
                 write(*,*)'         in calling program'
                 write(*,*)'         Number of visible triangles '
                 write(*,*)'         for this point             =',nvis
                 write(*,*)'         Current value of nv_max    =',nv_max
                 stop
              end if

	      do 60 i=1,nvis
                 t = vis_tlist(i)
                 ipos = vis_elist(i)
                 jpos = mod(vis_elist(i),3)+1
c	         write(6,*)' visible t =',t,' node',v(ipos,t),v(jpos,t)
 	         call addpoint
     &           (points,e,v,p,t,ipos,numtri,newpoint,add_tlist) 
                 newpoint = .false.
 60           continue

	   else
 
c	      write(6,*)' point located in triangle',t

c					add node to inside of convex hull
c					using swapping algorithm

  	      call insertpoint(points,e,v,p,t,numtri,iface) 

           end if

           if (numtri.gt.numtri_max) then
              write (*,*) 'Error in subroutine delaun:'
              write (*,*) 'Too many triangles'
              write(*,*)' Remedy: increase size of parameter numtri_max'
              write(*,*)'         in calling program'
              write(*,*)'         Number of triangles '
              write(*,*)'         for this point             =',numtri
              write(*,*)'         Current value of numtri_max    =',
     &                  numtrimax_max
              stop
           endif

 100    continue

	return
	end
