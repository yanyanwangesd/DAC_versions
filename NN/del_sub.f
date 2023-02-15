c------------------------------------------------------------------------
c
c	Subroutine visiblelist - calculates all sides of triangles 
c		                 visible from the point (x,y), which is
c			         outside of the convex hull.
c			
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		x,y			Co-ordinates of test point p
c		t			Index of any triangle on hull 
c					that is visible from point p.
c					(Usually given by routine Triloc_del.)
c		tpos			Position of edge in triangle t
c					(using Sloan's adjacency convention)
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c		nvis			Number of triangles visible from point p
c		vis_tlist		List of triangles visible from p
c		vis_elist		List of edges visible from p
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c		 This routine was converted from using Sloan's definition of
c		 the adjacency matrix to the `opposite' definition on 30/1/96.
c
c	Calls routine visible.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c		
	Subroutine visiblelist
     &             (points,neighbour,vertices,x,y,t,tpos,eps,
     &              vis_tlist,vis_elist,nvis)

	real*8		points(2,*)
	real*8		x,y
	real*8		eps
 	integer		vertices(3,*)
 	integer		neighbour(3,*)
 	integer		vis_tlist(*),vis_elist(*)
 	integer		t,tpos,pos,t1,t2,edg,tnew
	logical		visible
	logical		special
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	nvis = 1
	vis_tlist(1) = t
	vis_elist(1) = tpos
    	inode = vertices(tpos,t)
cd      write(6,100)t,inode,vertices(mod(tpos,3)+1,t)
    	pos = c1(tpos)
      	jnode = vertices(pos,t)
    	t1 = neighbour(c2(pos),t)
        special = .false.
	if(t1.eq.0)then
	  t1 = t
          tnew = 0
          special = .true.
	end if

  5     continue
        if(.not.special)then
           pos = edg(t1,jnode,vertices)
           tnew = neighbour(c2(pos),t1)
cd         write(6,*)' tnew =',tnew,' t1',t1,' jnode',jnode
        end if
        special = .false.
        if(tnew.eq.0)then
  6        continue
  	   if(visible(x,y,points,vertices,t1,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t1
	      vis_elist(nvis) = pos
cd            write(6,100)t1,jnode,vertices(mod(pos,3)+1,t1)
	   else
cd	      write(6,200)t1,jnode,vertices(mod(pos,3)+1,t1)
              go to 10
	   end if
           pos = c1(pos)
	   jnode = vertices(pos,t1)
    	   tnew = neighbour(c2(pos),t1)
	   if(tnew.eq.0) go to 6
           t1 = tnew
	   go to 5 
        else
cd	   write(6,300)t1,jnode,vertices(mod(pos,3)+1,t1)
	   t1 = tnew
	   go to 5 
	end if

  10	jnode = inode
    	pos = c2(tpos)
    	t2 = neighbour(c2(pos),t)
        special = .false.
	if(t2.eq.0)then
	  t2 = t
          tnew = 0
          special = .true.
	end if

  15    continue
        if(.not.special)then
           pos = c2(edg(t2,jnode,vertices))
           tnew = neighbour(c2(pos),t2)
        end if
        special = .false.
        if(tnew.eq.0)then
  16       continue
  	   if(visible(x,y,points,vertices,t2,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t2
	      vis_elist(nvis) = pos
cd            write(6,100)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
	   else
cd	      write(6,200)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
              go to 20
	   end if
	   jnode = vertices(pos,t2)
    	   pos = c2(pos)
    	   tnew = neighbour(c2(pos),t2)
	   if(tnew.eq.0)go to 16
	   t2 = tnew
	   go to 15
        else
cd	   write(6,300)t2,jnode,vertices(mod(pos,3)+1,t2)
	   t2 = tnew
	   go to 15
	end if

 20	continue
	      
 100    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is visible')
 200    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not visible')
 300    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not on convex hull')

	return
	end
c
c------------------------------------------------------------------------
c
c	Function visible - determines whether the triangle t is visible
c		           from the point p on edge tpos.
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        t			Triangle to be tested
c	        tpos			Edge to be tested in triangle t 
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        visible			Logical: = true if edge is visible  
c	                                         = false if edge is not visible
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c	Calls no other routines.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Function visible(x,y,points,vertices,t,tpos,eps)

	real*8		points(2,*)
	real*8		del1,del2
	real*8		x,y
 	integer		vertices(3,*)
 	integer		t,tpos
	logical		visible
	real*8		eps,del
	integer		c1(3)
        save            c1
	data		c1/2,3,1/

        j = c1(tpos)
c						test edge tpos in triangle t
        i1 = vertices(tpos,t)
        i2 = vertices(j,t)
        del1 = (points(2,i1)-y)*(points(1,i2)-x)
        del2 = (points(1,i1)-x)*(points(2,i2)-y)
	del = del1-del2
        if(del.gt.eps)then
           visible = .true.
	else
           visible = .false.
	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	addpoint - inserts a point into an existing delaunay triangulation 
c		   when point is outside of triangulation (but attached to
c		   triangle t) using the stacking procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle on convex hull visible 
c					from input point 
c		numtri			number of triangles in current
c					triangulation.
c		tpos			position of start node in triangle t
c		tri			list of triangles visible from point p
c		newpoint		logical = true if t is the first
c					triangle on the hull visible from p
c
c	Output:
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	When newpoint = false then there are multiple triangles
c	from the new point to the convex hull, and addpoint must
c	be called once for each attached triangle. In this case
c	the initialization of the adjacency list includes the
c	neighouring triangles already processed by addpoint, i.e.
c	those from point p to the hull.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine addpoint (points,e,v,p,t,tpos,numtri,newpoint,tri) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t,tpos,ccw
	logical		swap
	integer		tri(*)
	logical		newpoint
        save 		ip,lp
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	if(newpoint)then
	   ip = 0
	   lp = 0
        end if

c			Add new node to existing triangulation

c					create new triangle
        numtri = numtri + 1
        v1 = v(tpos,t)
        v2 = v(c1(tpos),t)
        if(ccw(points(1,v1),
     &         points(1,v2),
     &         points(1,p),k).eq.-1)then
               itemp = v1
               v1 = v2
               v2 = itemp
        end if
        v(1,numtri) = p
        v(2,numtri) = v1
        v(3,numtri) = v2

c					initialize adjacency list including
c					neighbouring triangles attached
c					from the point to the hull.

        e(c2(1),numtri) = 0
        e(c2(2),numtri) = t
        e(c2(3),numtri) = 0

c				
        if(.not.newpoint)then
           do 10 j=1,lp
              k = tri(j)
              if(v(2,k).eq.v1)then
c                write(6,*)' v1 match with node 2'
c                write(6,*)' current triangle',numtri,' new',k
c                write(6,*)' nodes:',v(1,k),v(2,k),v(3,k) 
c                write(6,*)' e mat:',e(c2(1),k),e(c2(2),k),e(c2(3),k) 
                 e(c2(1),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v1)then
                 e(c2(1),numtri) = k
                 e(c2(3),k) = numtri
              end if
              if(v(2,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(3),k) = numtri
              end if
 10        continue
        end if

c
c					initialize stack

 	call stackinit

c                                       update adjacency list
c                                       for triangle on old boundary
        e(c2(tpos),t) = numtri


c					add new triangle on stack

	call push(numtri)

c					loop while stack is not empty

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl) 
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = l
	   end if

	   if(b.ne.0)then
              call push(r)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = r
           end if
	   if(c.ne.0) e(edg(c,l,e),c)=r

	else

c					record triangle attached to p
	   ip = ip + 1
           tri(ip) = l

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
        call stackflush()

	lp = ip

c       write(6,*)' Number of triangles attached to last point',ip
c1	write(6,*)(tri(i),i=1,ip)
c	write(6,*)' triangles attached to last point on hull'
c       do 100 i=1,ip
c          it=tri(i)
c          do 101 k=1,3
c             l=mod(k,3)+1
c             if(e(k,it).eq.0)then
c                write(6,*)' t',it,' edge ',v(k,it),v(l,it)
c             end if
c 101      continue 
c 100   continue 
c

	return
	end
c
c------------------------------------------------------------------------
c
c	insertpoint - inserts a point into an existing delaunay triangulation 
c		      (when new point is inside triangle t) using the stacking 
c		      procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle containing input point
c		numtri			number of triangles in current
c					triangulation.
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Output:
c
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The new point is assumed to be inside the convex hull of the
c	existing triangulation.
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine insertpoint (points,e,v,p,t,numtri,iface) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t
	logical		swap
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/


c					add new node to existing triangulation
        if(iface.eq.0)then
	   a = e(c2(1),t)
	   b = e(c2(2),t)
	   c = e(c2(3),t)
	   v1 = v(1,t)
	   v2 = v(2,t)
	   v3 = v(3,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+2	
	   e(c2(2),t) = a	
	   e(c2(3),t) = numtri+1	
        

c					create new triangles
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v2
	   v(3,numtri)=v3
	   e(c2(1),numtri)=t
	   e(c2(2),numtri)=b
	   e(c2(3),numtri)=numtri+1
	   numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=numtri-1
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        else
           j = iface
           k = c1(j)
           i = c1(k)
	   a = e(c2(i),t)
	   b = e(c2(j),t)
	   c = e(c2(k),t)
	   v1 = v(i,t)
	   v2 = v(j,t)
	   v3 = v(k,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+1
	   e(c2(2),t) = a	
	   e(c2(3),t) = 0 	

c					create new triangle
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=0
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        end if
  
c
c					initialize stack

 	call stackinit

c					add new triangles on stack
c					and update adjacency list

	if(a.ne.0)call push(t)

	if(b.ne.0)then
          e(edg(b,t,e),b)=numtri-1
          call push(numtri-1)
        end if

	if(c.ne.0)then
          e(edg(c,t,e),c)=numtri
          call push(numtri)
        end if
c					loop while stack is not empty

        if(a.eq.0.and.b.eq.0.and.c.eq.0)go to 100

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl)
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   end if
	   if(b.ne.0) call push(r)
	   if(c.ne.0) e(edg(c,l,e),c)=r

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
 100    continue
        call stackflush()

	return
	end
c
c------------------------------------------------------------------------
c
c	Function edg - finds edge in triangle l which is adjacent 
c		       to triangle k.
c
c		       (From Sloan 1987)
c
c------------------------------------------------------------------------
c
	Function edg(l,k,e)
c
	integer		l,k,i,e(3,*),edg
c
	do 10 i=1,3
	   if(e(i,l).eq.k)then
              edg = i
              return
           end if
 10     continue

	write(*,*)' ***Error in function edg***'
	write(*,*)' ***Triangles not adjacent***'
	write(*,*)' triangle = ',l,' looking for triangle',k

	stop
	end
c
c------------------------------------------------------------------------
c
c	logical function swap - checks to see if point p lies 
c			        inside circumcircle about points p1,p2,p3
c				using the algorithm of Cline and Renka
c				(see Sloan 1987).
c
c------------------------------------------------------------------------
c
	Function swap(p1,p2,p3,p)

	logical		swap

	real*8		p(2),p1(2),p2(2),p3(2)
	real*8		x13,y13,x23,y23,x1p,y1p,x2p,y2p
	real*8		cosa,cosb,sina,sinb

	x13=p1(1)-p3(1)
	y13=p1(2)-p3(2)
	x23=p2(1)-p3(1)
	y23=p2(2)-p3(2)
	x1p=p1(1)-p(1)
	y1p=p1(2)-p(2)
	x2p=p2(1)-p(1)
	y2p=p2(2)-p(2)

	cosa = x13*x23 + y13*y23
	cosb = x2p*x1p + y1p*y2p

	if((cosa.ge.0.d0).and.(cosb.ge.0.d0))then
            swap = .false.
        else if((cosa.lt.0.d0).and.(cosb.lt.0.d0))then
            swap = .true.
        else
            sina=x13*y23-x23*y13
            sinb=x2p*y1p-x1p*y2p
	    if((sina*cosb+sinb*cosa).lt.0.d0)then
                swap = .true.
            else
                swap = .false.
            end if
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc_del - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node co-ordinates	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c	        loc			first guess of triangle containing
c					(x, y).
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c	        k			index of face through which the
c					algorithm last passed (used by
c					routine visbilelist if out = .true.)
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a visible triangle
c		 on the hull, out is set to .true., and k is set to the
c		 index of the face of triangle loc visible from the input point
c		 (used as a starting point by the routine visiblelist)
c
c		 This version also returns the parameter iface. 
c		 If iface .ne. 0 then the input point is on the face of 
c		 triangle t between nodes iface and mod(iface,3)+1 and 
c		 it is also on the convex hull.
c
c		 A point is assumed to be on the edge (or its extension)
c		 between two nodes if it is inside the triangle at a 
c		 distance >= eps.
c
c		 Can be extended to higher dimensions using a similar
c		 stepping mechanism but without angular test.
c
c	         This routine was converted from using Sloan's definition of
c	         the adjacency matrix to the `opposite' definition on 30/1/96.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine Triloc_del
     &                       (x,y,points,vertices,neighbour,loc,eps,
     &                        out,k,iface)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2,del
	real*8		eps
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/
        logical         new

	out = .false.
        new = .true.
        ic = 0

 10     continue
c					point is outside convex hull
        if(out)return
        iface = 0

        do 20 i=1,3
	   j = c1(i)
c	   k = c1(j)
c					definition of adjacency matrix

c					use Sloan's 
c					definition of adjacency matrix
	   k = i

           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
           del = del1-del2
 	   if(dabs(del).le.eps)then
              iface = i
	   else if(del.gt.0.d0)then
	      if(neighbour(c2(k),loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(c2(k),loc)
	      end if
              if(.not.new.and.loc.eq.loc1)then
                 write(*,100) 
                 write(*,*)' Current triangle:',loc, 
     &           ' last three:',loc1,loc2,loc3
                 write(*,*)' New point      x:',x,' y:',y
                 write(*,*)' Triangle ',loc,
     &           ' v:',(vertices(j,loc),j=1,3),
     &           ' n:',(neighbour(c2(j),loc),j=1,3)
c                write(*,*)' del',del,' del1',del1,' del2',del2
                 write(*,101) 
                 stop
              end if
              if(new)then
                ic = ic + 1
                if(ic.eq.3)new = .false.
              end if
              loc1 = loc2
              loc2 = loc3
              loc3 = loc
	      go to 10
	   end if
 20     continue
	
c						check if input point is
c						on the convex hull
c
c MODIFICATION MADE BY JEAN on 14/1/2004
c	  if(neighbour(c2(iface),loc).ne.0)iface = 0

        if (iface.ne.0) then
	  if(neighbour(c2(iface),loc).ne.0)iface = 0
        endif

c       if(iface.ne.0)then
c          j = mod(iface,3)+1
c          jj = vertices(iface,loc)
c          kk = vertices(j,loc)
c          write(*,*)' point on triangle between nodes ',
c    &               jj,' and',kk
c          write(*,*)' point is on the convex hull'
c       end if

 100    format(/'Error in subroutine Triloc_del:',//
     &  ' Infinite loop detected in walking triangle algorithm',/,
     &  ' Probably due to rounding error creating a flat triangle'/)

 101    format(/1x,'Remedy: '/
     &  ' Either increase size of parameter eps in calling routine '/
     &  ' or re-order input points by running program nn_hull '/)

	return
	end
c
c------------------------------------------------------------------------
c
c	Function ccw - used to test the orientation of three points
c
c		Input : 	points p1,p2 and p3 (vectors 2x1)
c				(e.g. p(1,2) = x co-ordinate of p2)
c
c		Output: 	ccw,I
c
c		ccw    k
c	  	 1     0	:The direction p1,p2,p3 is ccw (+ve)   
c	  	-1     0	:The direction p1,p2,p3 is  cw (-ve)   
c	  	 1     1	:p1,p2,p3 are colinear & p2 in middle  
c	  	-1     1	:p1,p2,p3 are colinear & p1 in middle
c	  	 0     1	:p1,p2,p3 are colinear & p3 in middle 
c
c
c				Calls no other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
      Integer Function ccw(p1,p2,p3,k)
c     
      real*8		p1(2),p2(2),p3(2)
      real*8		dx1,dx2,dy1,dy2,a,b
c      integer		ccw
c
      dx1 = p2(1) - p1(1)
      dx2 = p3(1) - p1(1)
      dy1 = p2(2) - p1(2)
      dy2 = p3(2) - p1(2)
      a = dx1*dy2
      b = dy1*dx2
      if (a.gt.b)then
         k = 0
         ccw = 1
      else if(a.lt.b)then
         k = 0
         ccw = -1
      else if(dx1*dx2.lt.0.0.or.dy1*dy2.lt.0.0)then
         k = 1
         ccw = -1 
      else if((dx1*dx1+dy1*dy1).lt.(dx2*dx2+dy2*dy2))then
         k = 1
         ccw = 1 
      else
         k = 1
         ccw = 0 
      end if
      return
      end
c
c------------------------------------------------------------------------
c
c	function theta - returns a real number between 0 and 360
c		         which has the same ordering as the angle
c			 between the line (a,b) and the horizontal.
c
c------------------------------------------------------------------------
c
	function theta(a,b)

	real*8		a(2),b(2)
	real*4		theta

	dx = b(1) - a(1)  
	ax = abs(dx)
	dy = b(2) - a(2)  
	ay = abs(dy)

	theta = 0.
	d = ax+ay
	if(d.ne.0)theta=dy/d

	if(dx.lt.0.0)then
           theta = 2.-theta
	else if(dy.lt.0.0)then
           theta = 4.+theta
	end if
	theta = theta*90.

	return
	end
