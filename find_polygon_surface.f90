subroutine find_polygon_surface (geometry,network,params,delaunay)
! Utility routine to compute the surface area associated with each node
! based on the real location of the divides around the node.
! More accurate and faster than the Voronoi cell calculation.
! Surface area is the sum of triangles surrounding a node.
! Each vertices of each triangle is the node it self, a divide, and
! a centroid of a delaunay tribagle.

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  integer i,j,k,h,tri, i1, i2, i3, h1, h2, tri1, tri2, div1, found
  double precision mid_x,mid_y,area, tri_area, total_area, total_area_2, part_area
  integer,dimension(:,:),allocatable::tri_per_node
  integer, dimension(:),allocatable::n_tri_per_node
  integer,dimension(:,:),allocatable::div_per_node
  integer, dimension(:),allocatable::n_div_per_node
  double precision, dimension(:),allocatable::area_per_triangle


  call time_in ('find_polygon_surface')
  allocate(tri_per_node(geometry%nnode,geometry%nnmax), n_tri_per_node(geometry%nnode))
  allocate(div_per_node(geometry%nnode,geometry%nnmax), n_div_per_node(geometry%nnode))
  allocate(area_per_triangle(delaunay%ntriangles))

  area_per_triangle = 0.d0
  tri_per_node = 0
  n_tri_per_node = 0
  do i = 1,delaunay%ntriangles
     i1 = delaunay%icon(1,i)
     i2 = delaunay%icon(2,i)
     i3 = delaunay%icon(3,i)
     n_tri_per_node(i1) =  n_tri_per_node(i1) + 1
     tri_per_node(i1,n_tri_per_node(i1)) = i
     n_tri_per_node(i2) =  n_tri_per_node(i2) + 1
     tri_per_node(i2,n_tri_per_node(i2)) = i
     n_tri_per_node(i3) =  n_tri_per_node(i3) + 1
     tri_per_node(i3,n_tri_per_node(i3)) = i
  enddo
  div_per_node = 0
  n_div_per_node = 0
  do i = 1,geometry%ndivide
     i1 = geometry%nndivnode(i,1)
     i2 = geometry%nndivnode(i,2)
     n_div_per_node(i1) =  n_div_per_node(i1) + 1
     div_per_node(i1,n_div_per_node(i1)) = i
     n_div_per_node(i2) =  n_div_per_node(i2) + 1
     div_per_node(i2,n_div_per_node(i2)) = i
  enddo





  total_area = 0.d0

  ! do some preprocessing so that the calculation is faster

  do i=1,geometry%nnode
     ! a node is in its polygon
     area = 0.d0
     do k = 1,geometry%nb(i)
        j = geometry%nn(k,i)
        if (i.eq.network%receiver(j)) then !i and j are part of the channel network
           ! find the centroids of the two triangles that share i and j 
           found = 0
           do h1 = 1,n_tri_per_node(i)
              tri1 = tri_per_node(i,h1)
              do h2 = 1,n_tri_per_node(j)
                 if (tri_per_node(j,h2).eq.tri1) then ! found one triangle
                    part_area = tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                         &delaunay%centers(1,tri1),delaunay%centers(2,tri1))
                    area = area + part_area
                    area_per_triangle(tri1) = area_per_triangle(tri1) + part_area
                    ! to find the other triangle go over the neighbors of tri1
                    tri2 = delaunay%neighbours(1,tri1)
!!$                     print*, 'coordinates of', tri1
!!$                     print*, '1', geometry%x(delaunay%icon(1,tri1)),geometry%y(delaunay%icon(1,tri1))
!!$                     print*, '2', geometry%x(delaunay%icon(2,tri1)),geometry%y(delaunay%icon(2,tri1))
!!$                     print*, '3', geometry%x(delaunay%icon(3,tri1)),geometry%y(delaunay%icon(3,tri1))
!!$                     print*, 'neighbors of',tri1
!!$                     print*, '1', delaunay%neighbours(1,tri1)
!!$                     print*, '2', delaunay%neighbours(2,tri1)
!!$                     print*, '3', delaunay%neighbours(3,tri1)
!!$                     print*, 'i and j are',i,j
                    if (tri2.ne.0) then
                       if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                          if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                             part_area = tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                  &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                             area = area + part_area
                             area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                             found = 1
                             exit
                          endif
                       endif
                    endif
                    tri2 = delaunay%neighbours(2,tri1)
                    if (tri2.ne.0) then
                       if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                          if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                             part_area = tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                  &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                             area = area + part_area
                             area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                             found = 1
                             exit
                          endif
                       endif
                    endif
                    tri2 = delaunay%neighbours(3,tri1)
                    if (tri2.ne.0) then
                       if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                          if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                             part_area = tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                  &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                             area = area + part_area
                             area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                             found = 1
                             exit
                          endif
                       endif
                    endif
                    print*, 'I should not be here - only one neighboring triangle'
                 endif
              enddo
              if (found.eq.1) exit
           enddo
        else
           !if i and j are not part of the channel network, must be divide in between or they are both at the base level
           if (j.ne.network%receiver(i)) then 
              if (geometry%fix(i).ge.1.and.geometry%fix(j).ge.1) then !both on the base level no divide in between
                 !find the common triangles
                 mid_x = (geometry%x(i) + geometry%x(j))/2.d0
                 mid_y = (geometry%y(i) + geometry%y(j))/2.d0
                 found = 0
                 do h1 = 1,n_tri_per_node(i)
                    tri1 = tri_per_node(i,h1)
                    do h2 = 1,n_tri_per_node(j)
                       if (tri_per_node(j,h2).eq.tri1) then ! found one triangle
                          part_area = tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                               &delaunay%centers(1,tri1),delaunay%centers(2,tri1))
                          area = area +  part_area
                          area_per_triangle(tri1) = area_per_triangle(tri1) + part_area
                          !looking for the next one
                          tri2 = delaunay%neighbours(1,tri1)
!!$                          print*, 'coordinates of', tri1
!!$                          print*, '1', geometry%x(delaunay%icon(1,tri1)),geometry%y(delaunay%icon(1,tri1))
!!$                          print*, '2', geometry%x(delaunay%icon(2,tri1)),geometry%y(delaunay%icon(2,tri1))
!!$                          print*, '3', geometry%x(delaunay%icon(3,tri1)),geometry%y(delaunay%icon(3,tri1))
!!$                          print*, 'neighbors of',tri1
!!$                          print*, '1', delaunay%neighbours(1,tri1)
!!$                          print*, '2', delaunay%neighbours(2,tri1)
!!$                          print*, '3', delaunay%neighbours(3,tri1)
!!$                          print*, 'i and j are',i,j
!!$                          print*, delaunay%icon(1,tri2),delaunay%icon(2,tri2),delaunay%icon(3,tri2)
                          if (tri2.ne.0) then
                             if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                                if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                                   part_area = tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                        &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                                   area = area + part_area
                                   area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                                   found = 1
                                   exit
                                endif
                             endif
                          endif
                          tri2 = delaunay%neighbours(2,tri1)
                          if (tri2.ne.0) then
                             if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                                if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                                   part_area = tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                        &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                                   area = area + part_area
                                   area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                                   found = 1
                                   exit
                                endif
                             endif
                          endif
                          tri2 = delaunay%neighbours(3,tri1)
                          if (tri2.ne.0) then
                             if (delaunay%icon(1,tri2).eq.i.or.delaunay%icon(2,tri2).eq.i.or.delaunay%icon(3,tri2).eq.i) then
                                if (delaunay%icon(1,tri2).eq.j.or.delaunay%icon(2,tri2).eq.j.or.delaunay%icon(3,tri2).eq.j) then
                                   part_area = tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                        &delaunay%centers(1,tri2),delaunay%centers(2,tri2))
                                   area = area + part_area
                                   area_per_triangle(tri2) = area_per_triangle(tri2) + part_area
                                   found = 1
                                   exit
                                endif
                             endif
                          endif
                          exit ! it is possible that there is only one common triangle
                       endif
                    enddo
                    if (found.eq.1) exit
                 enddo
              else !must be a divide between them, find it
                 found = 0
                 do h1 = 1,n_div_per_node(i)
                    div1 = div_per_node(i,h1)
                    do h2 = 1,n_div_per_node(j)
                       if (div_per_node(j,h2).eq.div1) then
                          !print*, i,j,div1, geometry%nndivtri(div1,1), geometry%nndivtri(div1,2)
                          !print*, n_div_per_node(i),n_div_per_node(j)
                          tri = geometry%nndivtri(div1,1)
                          part_area = tri_area(geometry%x(i),geometry%y(i),geometry%xdiv(div1),geometry%ydiv(div1),&
                               &delaunay%centers(1,tri),delaunay%centers(2,tri)) 
                          area = area + part_area
                          area_per_triangle(tri) = area_per_triangle(tri) + part_area

                          tri = geometry%nndivtri(div1,2)
                          part_area = tri_area(geometry%x(i),geometry%y(i),geometry%xdiv(div1),geometry%ydiv(div1),&
                               &delaunay%centers(1,tri),delaunay%centers(2,tri)) 
                          area = area +  part_area
                          area_per_triangle(tri) = area_per_triangle(tri) + part_area
                          found = 1
                          exit !because the divide was found
                       endif
                    enddo
                    if (found.eq.1) exit
                 enddo
              endif
           endif
        endif
     enddo
     geometry%surface(i) = area
     total_area = total_area + area
     !print*, 'area is', area
  enddo


!!$  print*, 'total area is', total_area
!!$  total_area_2 = 0.
!!$  do i=1,delaunay%ntriangles
!!$     i1 = delaunay%icon(1,i)
!!$     i2 = delaunay%icon(2,i)
!!$     i3 = delaunay%icon(3,i)
!!$     part_area = tri_area(geometry%x(i1),geometry%y(i1),geometry%x(i2),geometry%y(i2),geometry%x(i3),geometry%y(i3))
!!$     if (abs((part_area-area_per_triangle(i))/part_area).gt.0.01)then
!!$        print*,'area problem for triangle',i
!!$        print*, 'real area=',part_area,'calculated area=',area_per_triangle(i)
!!$        print*, 'coordinates of',i
!!$        print*, '1', geometry%x(delaunay%icon(1,i)),geometry%y(delaunay%icon(1,i)),geometry%z(delaunay%icon(1,i))
!!$        print*, '2', geometry%x(delaunay%icon(2,i)),geometry%y(delaunay%icon(2,i)),geometry%z(delaunay%icon(2,i))
!!$        print*, '3', geometry%x(delaunay%icon(3,i)),geometry%y(delaunay%icon(3,i)),geometry%z(delaunay%icon(3,i))
!!$     endif
!!$     total_area_2 = total_area_2 + part_area
!!$  enddo
!!$  print*, 'area according to triangulation is', total_area_2

  deallocate(tri_per_node,n_tri_per_node,div_per_node,n_div_per_node,area_per_triangle)
  call time_out ('find_polygon_surface')
  return

end subroutine find_polygon_surface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function tri_area (x1,y1,x2,y2,x3,y3)
  implicit none
  double precision x1,y1,x2,y2,x3,y3,l1,l2,l3,s
  
  l1 = dsqrt((x1-x2)**2.d0 + (y1-y2)**2.d0)
  l2 = dsqrt((x2-x3)**2.d0 + (y2-y3)**2.d0)
  l3 = dsqrt((x3-x1)**2.d0 + (y3-y1)**2.d0)
 ! s = (l1+l2+l3)/2.
  !tri_area =sqrt(s*(s-l1)*(s-l2)*(s-l3))  
  tri_area = 0.5d0*abs(x1*(y3-y2)+x2*(y1-y3) + x3*(y2-y1))
  
  return 
end function tri_area
