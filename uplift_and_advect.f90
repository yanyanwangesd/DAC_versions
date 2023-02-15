subroutine uplift_and_advect (geometry,params)

  ! Update the node geometry (for both horizontal and vertical tectonic advection)
  ! Note that horizontal advection may be turned off by setting the params%move_points flag to .FALSE.
  ! in the initialize_parameters routine

  use definitions

  implicit none

  type (geom) geometry
  type (parm) params
  integer i,j,k, ncount, bcon, corner, nbcon
  double precision y1,y2,y3,y4,x1,x2, nelev
  double precision xmin,xmax,ymin,ymax, xsc
  double precision,allocatable:: bold(:)

  call time_in ('uplift_and_advect')

  !when advecting nodes there are special add and remove operations for which the following is needed:
  if (params%move_points .or. 1.eq.1) then
     xmin = minval(geometry%x)
     xmax = maxval(geometry%x)
     ymin = minval(geometry%y)
     ymax = maxval(geometry%y)
  endif



  do i=1,geometry%nnode
     ! non-boundary nodes
     if(geometry%fix(i).eq.0 ) then
        geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
         if(params%move_points) then
            geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i)
            geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
         endif
     else
     !boundary nodes
     !0     z=0
     !1     z=zinit
     !2     z=zneigh
     !3     -
     !4     z=z+w*dt
     !5     -
     !6     -
     !7     -
     !8     x=x+u*dt, z=0
     !9     x=x+u*dt, z=zinit
     !10    x=x+u*dt, z=neigh
     !11    -
     !12    x=x+u*dt, z=z+w*dt
     !13    -
     !14    -
     !15    -
     !16    y=y+v*dt, z=0
     !17    y=y+v*dt, z=zinit
     !18    y=y+v*dt, z=zneigh
     !19    -
     !20    y=y+v*dt, z=z+w*dt
     !21    -
     !22    -
     !23    -
     !24    x=x+u*dt, y=y+v*dt, z=0
     !25    x=x+u*dt, y=y+v*dt, z=zinit
     !26    x=x+u*dt, y=y+v*dt, z=zneigh
     !27    -
     !28    x=x+u*dt, y=y+v*dt, z=z=z+w*dt
       ! LEFT BOUNDARY
         if(geometry%fix(i).eq.1)then
            ! z
            select case ( abs(geometry%boundary(1))  )
                case(0,8,16,24)
                    geometry%z(i)=0.d0
                case(1,9,17,25)
                     !geometry%z(i)=geometry%z(i)
                !case(2,10,18,26)
                case(4,12,20,28)
                     geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
            end select
            if(params%move_points) then
                select case ( abs(geometry%boundary(1))  )
                    case(8,9,10,12,24, 25,26,28)
                        geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i)
                end select
                select case ( abs(geometry%boundary(1))  )
                    case(16,17,18,20,24,25,26,28)
                         geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
                end select
            endif
         endif
         !RIGHT BOUNDARY
          if(geometry%fix(i).eq.2)then
            ! z
            select case (  abs(geometry%boundary(2))  )
                case(0,8,16,24)
                    geometry%z(i)=0.d0
                case(1,9,17,25)
                     !geometry%z(i)=geometry%z(i)
                !case(2,10,18,26)
                case(4,12,20,28)
                     geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
            end select
            if(params%move_points) then
                select case (  abs(geometry%boundary(2))  )
                    case(8,9,10,12,24, 25,26,28)
                        geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i)
                end select
                select case (  abs(geometry%boundary(2))  )
                    case(16,17,18,20,24,25,26,28)
                         geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
                end select
            endif
         endif
         ! FRONT BOUNDARY
         if(geometry%fix(i).eq.3)then
            ! z
            select case ( abs(geometry%boundary(3))  )
                case(0,8,16,24)
                    geometry%z(i)=0.d0
                case(1,9,17,25)
                     !geometry%z(i)=geometry%z(i)
                !case(2,10,18,26)
                case(4,12,20,28)
                     geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
            end select
            if(params%move_points) then
                select case ( abs(geometry%boundary(3))  )
                    case(8,9,10,12,24, 25,26,28)
                        geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i)
                end select
                select case (  abs(geometry%boundary(3))  )
                    case(16,17,18,20,24,25,26,28)
                         geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
                end select
            endif
         endif
         ! BACK BOUNDARY
                  if(geometry%fix(i).eq.4)then
            ! z
            select case (abs(geometry%boundary(4))  )
                case(0,8,16,24)
                    geometry%z(i)=0.d0
                case(1,9,17,25)
                     !geometry%z(i)=geometry%z(i)
                !case(2,10,18,26)
                case(4,12,20,28)
                     geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
            end select
            if(params%move_points) then
                select case ( abs(geometry%boundary(4))  )
                    case(8,9,10,12,24, 25,26,28)
                        geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i)
                end select
                select case ( abs(geometry%boundary(4) ) )
                    case(16,17,18,20,24,25,26,28)
                         geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
                end select
            endif
         endif

     endif
  enddo





!! ---- set boundary nodes to averages of neighbours
 do i=1,geometry%nnode
    if(geometry%fix(i).ge.1 )then
         bcon = abs(geometry%boundary( geometry%fix(i) ))
         if(  bcon.eq.2 .or. bcon.eq.10.  .or. bcon.eq.18  .or. bcon.eq.26   ) then
             ncount=0
             nelev=0
             !find neighbours that are not boundaries
             do k=1,geometry%nb(i)
                j=geometry%nn(k,i)
                if( geometry%fix(j).eq.0 )then
                        ncount=ncount+1
                        nelev=nelev+geometry%z(j)
                endif
             enddo
             ! in case no non-boundary neighbours are found - sometimes corners -  set it from all neighbours
             if (ncount.eq.0)then
                 do k=1,geometry%nb(i)
                    j=geometry%nn(k,i)
                    ncount=ncount+1
                    nelev=nelev+geometry%z(j)
                 enddo
             endif
             !do average
             if(ncount.gt.0) geometry%z(i)=nelev/ncount
             ! do adjust only to part of the difference
             !if(ncount.gt.0) geometry%z(i)=    geometry%z(i) *0.33d0 + 0.67d0* nelev/ncount  - 0.01d0
         endif
    endif
 enddo
 !------





  call time_out ('uplift_and_advect')

  return

end subroutine uplift_and_advect
