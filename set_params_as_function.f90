subroutine set_params_as_function (geometry,params)

  ! Routine to ...

  use definitions

  implicit none

  type (geom) geometry
  type (parm) params

  integer i,j,k,l, nset, nsetr, icond, inbound, ri,rj
  double precision  d1,d2,d3, fmin,fmax, fwk, foffs, tmin,tmax, rastx,rasty, rxsize, rysize, rx0, ry0, wsum, valsum, dist,distx,disty
  double precision,pointer:: fdep(:), depvar(:)

!   ! legacy line of code that is needed so that find_precipitation does not need to be called  - set anyway where dependent variable is not P
!   geometry%precipitation=params%rainfall_height
!   ! and the resetting part from uplift_and_advect
!   geometry%u = 0.
!   geometry%v = 0.


  if(.NOT.params%f_varies) return

  do nset=1, params%f_num_sets

      ! dependent variable, select, and set pointer to
      select case(params%f_variable_determined(nset))
        case(0)
            depvar=>geometry%u
        case(1)
            depvar=>geometry%v
        case(2)
            depvar=>geometry%w
        case(3)
            depvar=>geometry%precipitation
        case(4)
            depvar=>geometry%k
        case(5)
            depvar=>geometry%z
      end select

      ! polynome or raster
      if(params%f_ctype(nset).eq.1)then
      !!!!!! this is a polynomial condition !!!!!!!!!!!!!!!!!!!!!!!!!!
            ! independent variable, select, and set pointer to
            select case(params%f_depends_on(nset))
              case(0)
                  fdep=>geometry%x
              case(1)
                  fdep=>geometry%y
              case(2)
                  fdep=>geometry%z
            end select

              tmin = params%f_timebound(nset,1)
              tmax = params%f_timebound(nset,2)
              if(tmin.eq.0 .and. tmax.eq.0) tmax = 1.d18

              if(  params%time.ge.tmin  .and. params%time.le.tmax  )then
                      do i = 1, geometry%nnode
                          ! which condition segment (line) applies
                          inbound=0
                          do icond=1,params%f_polyseg(nset)
                              fmin=params%poly(nset,icond,1)
                              fmax=params%poly(nset,icond,2)
                              foffs=params%poly(nset,icond,3)
                              if(fdep(i).ge.fmin  .and. fdep(i).le.fmax)then
                                  inbound=1
                                  exit
                              endif
                          enddo
                          if(inbound.eq.0)then                            !!RECONSIDER - simply not setting when out of bounds
                              if(fdep(i).lt.params%poly(nset,1,1))then
                                    icond=1
                                    !print*, "node below lower bound of condition coverage"
                              else if(fdep(i).gt.params%poly(nset,params%f_polyseg(nset),2))then
                                    icond=params%f_polyseg(nset)
                                    !print*, "node above upper bound of condition coverage"
                              else
                                    stop "fatally screwed up polynome input - this can be the case if defined segments do not connect"
                              endif
                          endif
                          fwk=fdep(i)-foffs

                          d3=0
                          do j=1,params%pn(nset,icond)  ! loop and sum up contributions of all polynomial degrees
                              d3=d3  + params%poly(nset,icond,j+3) * (fwk)**(j-1)
                          enddo
                            ! this is a boundary, and we set topography
                            if(geometry%fix(i).ne.0   .and. params%f_variable_determined(nset).eq.5)then
                                ! fixed to 0 m - 0,8,16,24
                                if( mod ( geometry%boundary(geometry%fix(i)) ,8)  .eq. 0  )then
                                      geometry%z(i)=0
                                ! fixed to initial height - 1,9,17,25
!                                else if( mod ( geometry%boundary(geometry%fix(i)) ,8)  .eq. 1  )then
!                                      geometry%z(i)=d3
                                ! types of boundaries that move vertically
                                else
                                      select case(params%f_superpose(nset))
                                          case(0)
                                              depvar(i)=d3
                                          case(1)
                                              depvar(i)=depvar(i)+d3
                                          case(2)
                                              depvar(i)=depvar(i)*d3
                                      end select
                               endif
                            ! internal node, or we don't set topography
                            else
                                  select case(params%f_superpose(nset))
                                      case(0)
                                          depvar(i)=d3
                                      case(1)
                                          depvar(i)=depvar(i)+d3
                                      case(2)
                                          depvar(i)=depvar(i)*d3
                                  end select
                            endif
                      enddo
              endif

              ! dissociate pointers for next condition
              nullify(fdep)
              nullify(depvar)
          !!!!!! this is a polynomial condition !!!!!!!!!!!!!!!!!!!!!!!!!!


        else
          !!!!!! this is a raster condition !!!!!!!!!!!!!!!!!!!!!!!!!!

              ! recover raster local index from global condition index
              nsetr = params%f_cidx(nset,2)

              tmin = params%f_timebound(nset,1)
              tmax = params%f_timebound(nset,2)
              if(tmin.eq.0 .and. tmax.eq.0) tmax = 1.d18

              if(  params%time.ge.tmin  .and. params%time.le.tmax  )then
                  !determine extent and resolution of grid
                  rxsize = geometry%xl + params%f_rmarg(nsetr,1) + params%f_rmarg(nsetr,2)
                  rysize = geometry%yl + params%f_rmarg(nsetr,3) + params%f_rmarg(nsetr,4)
                  rx0 = -params%f_rmarg(nsetr,1)
                  ry0 = -params%f_rmarg(nsetr,3)
                  rastx = rxsize/( params%f_rnnode(nsetr,2)  -1)
                  rasty = rysize/( params%f_rnnode(nsetr,1)  -1)

                  rastx = rxsize/( params%f_rnnode(nsetr,2)  )
                  rasty = rysize/( params%f_rnnode(nsetr,1)  )

                  do i = 1, geometry%nnode

                    ! find closest raster point - dist with reference to (shifted) origin
                    ! nearest is cheaper, interpolation from floor more accurate
                    rj = floor ( ( geometry%x(i) + params%f_rmarg(nsetr,1) )/rastx ) +1
                    ri = floor ( ( geometry%y(i) + params%f_rmarg(nsetr,3) )/rasty ) +1
                    if(rj.lt.1)rj = 1
                    if(rj.gt.params%f_rnnode(nsetr,2))rj = params%f_rnnode(nsetr,2)-1
                    if(ri.lt.1)ri = 1
                    if(ri.gt.params%f_rnnode(nsetr,1))ri = params%f_rnnode(nsetr,1)-1

                    d3 = params%f_raster(nsetr)%matx(ri,rj)


!                     ! alternatively do interpolation
!                     rj = floor ( ( geometry%x(i) + params%f_rmarg(nsetr,1) )/rastx ) +1
!                     ri = floor ( ( geometry%y(i) + params%f_rmarg(nsetr,3) )/rasty ) +1
!                     if(rj.lt.1)rj = 1
!                     if(rj.ge.params%f_rnnode(nsetr,2))rj = params%f_rnnode(nsetr,2)-1
!                     if(ri.lt.1)ri = 1
!                     if(ri.ge.params%f_rnnode(nsetr,1))ri = params%f_rnnode(nsetr,1)-1

!                     wsum = 0
!                     valsum = 0

!                     distx =  geometry%x(i) - ( (rj-1)*rastx - params%f_rmarg(nsetr,1) )
!                     disty =  geometry%y(i) - ( (ri-1)*rasty - params%f_rmarg(nsetr,3) )
!                     dist = sqrt(  distx**2 + disty**2   )
!                     !dist = (distx) * (disty)
!                     wsum  = wsum + 1/dist
!                     valsum = valsum + params%f_raster(nsetr)%matx(ri,rj) * 1/dist

!                     distx =  geometry%x(i) - ( (rj)*rastx - params%f_rmarg(nsetr,1) )
!                     disty =  geometry%y(i) - ( (ri-1)*rasty - params%f_rmarg(nsetr,3) )
!                     dist = sqrt(  (distx)**2 + disty**2   )
!                     !dist = (rastx-distx) * (disty)
!                     wsum  = wsum + 1/dist
!                     valsum = valsum + params%f_raster(nsetr)%matx(ri,rj+1) * 1/dist

!                     distx =  geometry%x(i) - ( (rj-1)*rastx - params%f_rmarg(nsetr,1) )
!                     disty =  geometry%y(i) - ( (ri)*rasty - params%f_rmarg(nsetr,3) )
!                     dist = sqrt(  distx**2 + (disty)**2   )
!                     !dist = (distx) * (rasty-disty)
!                     wsum  = wsum + 1/dist
!                     valsum = valsum + params%f_raster(nsetr)%matx(ri+1,rj) * 1/dist

!                     distx =  geometry%x(i) - ( (rj)*rastx - params%f_rmarg(nsetr,1) )
!                     disty =  geometry%y(i) - ( (ri)*rasty - params%f_rmarg(nsetr,3) )
!                     dist = sqrt(  (distx)**2 + (disty)**2   )
!                     !dist = (rastx - distx) * (rasty - disty)
!                     wsum  = wsum + 1/dist
!                     valsum = valsum + params%f_raster(nsetr)%matx(ri+1,rj+1) * 1/dist
!                     d3 = valsum / wsum
!                    !


                    ! this is a boundary, and we set topography
                    if(geometry%fix(i).ne.0   .and. params%f_variable_determined(nset).eq.5)then
                        ! fixed to 0 m - 0,8,16,24
                        if( mod ( geometry%boundary(geometry%fix(i)) ,8)  .eq. 0  )then
                              geometry%z(i)=0
                        ! fixed to initial height - 1,9,17,25
!                        else if( mod ( geometry%boundary(geometry%fix(i)) ,8)  .eq. 1  )then
!                              geometry%z(i)=d3
                        ! types of boundaries that move vertically
                        else
                              select case(params%f_superpose(nset))
                                  case(0)
                                      depvar(i)=d3
                                  case(1)
                                      depvar(i)=depvar(i)+d3
                                  case(2)
                                      depvar(i)=depvar(i)*d3
                              end select
                       endif
                    ! internal node, or we don't set topography
                    else
                          select case(params%f_superpose(nset))
                              case(0)
                                  depvar(i)=d3
                              case(1)
                                  depvar(i)=depvar(i)+d3
                              case(2)
                                  depvar(i)=depvar(i)*d3
                          end select
                    endif

                  enddo  ! node loop
              endif
              ! dissociate pointers for next condition
              nullify(depvar)
          !!!!!! this is a raster condition !!!!!!!!!!!!!!!!!!!!!!!!!!

    endif



  enddo   ! loop all sets

end subroutine set_params_as_function
