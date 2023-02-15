subroutine reboundary(geometry,params)


  use definitions

  implicit none

  type (parm) params
  type (geom) geometry


  integer i,j, k, ncount, bcon, corner, nbcon,ccount
  double precision y1,y2,y3,y4,x1,x2, nelev
  double precision xmin,xmax,ymin,ymax, xsc, xarr(3), yarr(3)
  double precision,allocatable:: bold(:)


! extra provision to make sure that corners behave as desired
  ! find limiting boundary condition, and toggle the corner's boundary affiliation
  !if(params%istep.le.5)then
    allocate(bold(geometry%nnode))
    bold=geometry%fix(1:geometry%nnode)
    ccount=0
    do i=1,geometry%nnode
        ! identify as boundary node
        if(bold(i).gt.0)then
            !try to find corner node
            bcon=bold(i)
            corner=0
            xarr=0
            yarr=0
            ncount=0
             do j=1, geometry%nb(i)
                k=geometry%nn(j,i)
                if(geometry%fix(k).gt.0)then
                    ncount=ncount+1
                    xarr(ncount)=geometry%x(k)
                    yarr(ncount)=geometry%y(k)
                endif
             enddo
            xmin=minval( xarr(1:ncount) )
            xmax=maxval(  xarr(1:ncount)  )
            ymin=minval(  yarr(1:ncount) )
            ymax=maxval(  yarr(1:ncount) )
            do j=1, geometry%nb(i)
                k=geometry%nn(j,i)
                if( bold(k).gt.0    .and. bold(k).eq.geometry%fix(k)  .and. bold(k).ne.bcon  .and. geometry%nb(i).lt.4  .and. (  xmax-xmin.gt.params%xc  .and.  ymax-ymin.gt.params%xc )  )then
                    corner=1
                    nbcon=bold(k)
                    exit
                endif
            enddo
            if(corner.eq.1)then
                !geometry%z(i)=9 !stupid marking
                select case(bcon)
                    !left
                    case(1)
                        ! left-back corner
                        if(geometry%y(i).gt.(ymax-ymin)/2)then
                            select case(geometry%boundary(4))
                                ! moving of back boundary in x allowed
                                !case(8,9,10,12,24, 25,26,28)
                                case(16,17,18,20,24, 25,26,28)
                                    !if(geometry%boundary(1).lt.8 .or.  (geometry%boundary(1).gt.12  .and. geometry%boundary(1).lt.24)  ) bcon=bcon!no problem
                                    if(geometry%boundary(1).lt.16) bcon=bcon!no problem
                                ! back boundary not moving in x
                                case default
                                    ! back boundary does not move, but left one does
                                    !if( (geometry%boundary(1).ge.8 .and. geometry%boundary(1).le.12)     .or.  geometry%boundary(1).ge.24) then
                                    if( geometry%boundary(1).ge.16 )then
                                        geometry%fix(i)=4    ! make the corner rigid/immobile
                                         bold(i)=-1
                                         ccount=ccount+1
                                    endif
                            end select
                        ! left-front corner
                        else
                            select case(geometry%boundary(3))
                                ! moving of front boundary in x allowed
                                !case(8,9,10,12,24, 25,26,28)
                                case(16,17,18,20,24, 25,26,28)
                                    !if(geometry%boundary(1).lt.8 .or.   (geometry%boundary(1).gt.12  .and. geometry%boundary(1).lt.24) ) bcon=bcon!no problem
                                    if(geometry%boundary(1).lt.16) bcon=bcon!no problem
                                ! front boundary not moving in x
                                case default
                                    ! front boundary does not move, but left one does
                                    !if( (geometry%boundary(1).ge.8 .and. geometry%boundary(1).le.12)     .or.  geometry%boundary(1).ge.24)then
                                    if( geometry%boundary(1).ge.16 )then
                                         geometry%fix(i)=3    ! make the corner rigid/immobile
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        endif
                    !right
                    case(2)
                        ! right-back corner
                        if(geometry%y(i).gt.(ymax-ymin)/2)then
                            select case(geometry%boundary(4))
                                ! moving of back boundary in x allowed
                                !case(8,9,10,12,24, 25,26,28)
                                case(16,17,18,20,24, 25,26,28)
                                    !if(geometry%boundary(2).lt.8 .or.   (geometry%boundary(1).gt.12  .and. geometry%boundary(1).lt.24) ) bcon=bcon!no problem
                                    if(geometry%boundary(2).lt.16) bcon=bcon!no problem
                                ! back boundary not moving in x
                                case default
                                    ! back boundary does not move, but right one does
                                    !if( (geometry%boundary(2).ge.8 .and. geometry%boundary(2).le.12)     .or.  geometry%boundary(2).ge.24)then
                                    if( geometry%boundary(2).ge.16 )then
                                         geometry%fix(i)=4    ! make the corner rigid/immobile
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        ! right-front corner
                        else
                            select case(geometry%boundary(3))
                                ! moving of front boundary in x allowed
                                !case(8,9,10,12,24, 25,26,28)
                                case(16,17,18,20,24, 25,26,28)
                                    !if(geometry%boundary(2).lt.8 .or.   (geometry%boundary(1).gt.12  .and. geometry%boundary(1).lt.24) ) bcon=bcon!no problem
                                    if(geometry%boundary(2).lt.16) bcon=bcon!no problem
                                ! front boundary not moving in x
                                case default
                                    ! front boundary does not move, but right one does
                                    !if( (geometry%boundary(2).ge.8 .and. geometry%boundary(2).le.12)     .or.  geometry%boundary(2).ge.24)then
                                    if( geometry%boundary(2).ge.16 )then
                                         geometry%fix(i)=3    ! make the corner rigid/immobile
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        endif
                    !front
                    case(3)
                        ! front-left corner
                        if(geometry%x(i).lt.(xmax-xmin)/2)then
                            select case(geometry%boundary(1))
                                ! moving of left boundary in y allowed
                                !case(16,17,18,20,24, 25,26,28)
                                case(8,9,10,12,24, 25,26,28)
                                    !if(geometry%boundary(3).lt.16) bcon=bcon!no problem
                                    if(geometry%boundary(3).lt.8 .or.  (geometry%boundary(3).gt.12  .and. geometry%boundary(3).lt.24)  ) bcon=bcon!no problem
                                ! left boundary not moving in y
                                case default
                                    ! left boundary does not move, but front one does
                                    !if( geometry%boundary(3).ge.16 )then
                                    if( (geometry%boundary(3).ge.8 .and. geometry%boundary(3).le.12)     .or.  geometry%boundary(3).ge.24) then
                                         geometry%fix(i)=1    ! make the corner rigid/
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        ! front-right corner
                        else
                            select case(geometry%boundary(2))
                                ! moving of right boundary in y allowed
                                !case(16,17,18,20,24, 25,26,28)
                                case(8,9,10,12,24, 25,26,28)
                                    !if(geometry%boundary(3).lt.16 ) bcon=bcon!no problem
                                    if(geometry%boundary(3).lt.8 .or.  (geometry%boundary(3).gt.12  .and. geometry%boundary(3).lt.24)  ) bcon=bcon!no problem
                                ! right boundary not moving in y
                                case default
                                    ! front boundary does not move, but left one does
                                    !if( geometry%boundary(3).ge.16)then
                                    if( (geometry%boundary(3).ge.8 .and. geometry%boundary(3).le.12)     .or.  geometry%boundary(3).ge.24) then
                                         geometry%fix(i)=2    ! make the corner rigid/immobile
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        endif
                    !back
                    case(4)
                        ! back-left corner
                        if(geometry%x(i).lt.(xmax-xmin)/2)then
                            select case(geometry%boundary(1))
                                ! moving of left boundary in y allowed
                                !case(16,17,18,20,24, 25,26,28)
                                case(8,9,10,12,24, 25,26,28)
                                    !if(geometry%boundary(4).lt.16) bcon=bcon!no problem
                                    if(geometry%boundary(4).lt.8 .or.  (geometry%boundary(4).gt.12  .and. geometry%boundary(4).lt.24)  ) bcon=bcon!no problem
                                ! left boundary not moving in y
                                case default
                                    ! left boundary does not move, but front one does
                                    !if( geometry%boundary(4).ge.16 )then
                                    if( (geometry%boundary(4).ge.8 .and. geometry%boundary(4).le.12)     .or.  geometry%boundary(4).ge.24) then
                                         geometry%fix(i)=1    ! make the corner rigid/immobile
                                          bold(i)=-1
                                          ccount=ccount+1
                                    endif
                            end select
                        ! back-right corner
                        else
                            select case(geometry%boundary(2))
                                ! moving of right boundary in y allowed
                                !case(16,17,18,20,24, 25,26,28)
                                case(8,9,10,12,24, 25,26,28)
                                   ! if(geometry%boundary(4).lt.16 ) bcon=bcon!no problem
                                    if(geometry%boundary(4).lt.8 .or.  (geometry%boundary(4).gt.12  .and. geometry%boundary(4).lt.24)  ) bcon=bcon!no problem
                                ! right boundary not moving in y
                                case default
                                    ! front boundary does not move, but left one does
                                    !if( geometry%boundary(4).ge.16)then
                                    if( (geometry%boundary(4).ge.8 .and. geometry%boundary(4).le.12)     .or.  geometry%boundary(4).ge.24) then
                                         geometry%fix(i)=2    ! make the corner rigid/immobile
                                         bold(i)=-1
                                         ccount=ccount+1
                                    endif
                            end select
                            !geometry%fix(i)=2
                        endif
                end select
            endif
        endif
    enddo
    deallocate(bold)
 ! endif

end subroutine reboundary
