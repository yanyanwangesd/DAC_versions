subroutine erode (geometry,network,stack,delaunay,params)

  ! Subroutine to compute fluvial erosion using a fully implicit method
  ! for n=1 or n=2
  ! Other cases not considered yet


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (parm) params
  type (del) delaunay
  integer i,is,j,lakeflag,icount, jcount,dep_nnode, k, l, m, icond
  integer next_node_ind,n,flag_not_in_list, index_of_receiver
  double precision fact,deltal,factp,erosion,sedheight,zi,zj,areai,areaj,length,depgrad,delz, erosion2, slopek
  double precision rtsafe, sedV, next_node_z, sumS, sum1, sum2, sum3, sum4, temp_e_rate,minsedV
  integer, dimension(:),allocatable::dep_index
  double precision, dimension(:),allocatable::ztest, areas, lengths, elevs,dz


  call time_in ('erode')

  ! set a depositional gradient for basins

  depgrad=.1d-6
  minsedV=.1d-4
  

  ! check for rivers flowing uphill
  do i = 1, geometry%nnode
     if(network%receiver(i).ne.0.and.geometry%fix(i).eq.0)then
        if(geometry%z(i).lt.geometry%z(network%receiver(i)))then
           print*, 'Warning River flowing uphill in erode.f90',geometry%z(i),geometry%z(network%receiver(i))
        endif
     endif
  enddo



  !NEW ALGORITHM FOR DEPOSITION IN LAKES. L.G.
  !check if the volume of sediments available for sedimentation 
  !can be deposited within a group of the lake nodes while maintaining depgrad slope
  !without reversing flow direction.
  !if yes - sediment it
  !if no - add another another node to the group which is the lowest in elevation donor of 
  ! any of the nodes of the group.   


  do i=1,geometry%nnode
     if(network%receiver(i).eq.0.and.geometry%fix(i).eq.0)then   !find lake nodes
        sedV=params%deltat*geometry%sediment_flux(i) !volume of sediments to be deposited in that lake
        if (network%ndon(i).eq.0) cycle
        if (sedV.eq.0..and.network%ndon(i).ne.0) then           
           sedV=minsedV
        endif
        
        !allocate data structures
        allocate(dep_index(geometry%nnode),ztest(geometry%nnode), areas(geometry%nnode), lengths(geometry%nnode), elevs(geometry%nnode),dz(geometry%nnode))
        dep_index=-1
        ztest=0.d0
        areas=0.d0
        lengths=0.d0
        elevs=0.d0
        dz=0.d0
        !first try to deposit within the node itself        
        dep_nnode=1
        dep_index(dep_nnode)=i
        ztest(dep_nnode)=geometry%z(i) + sedV/geometry%surface(i)
        areas(dep_nnode)=geometry%surface(i)
        lengths(dep_nnode)=0.d0
        elevs(dep_nnode)=geometry%z(i)

152     do j=1,dep_nnode !loop over the nodes that are being sedimenting
           k=dep_index(j) 
           if (network%ndon(k).ne.0) then
              do l=1,network%ndon(k) !loop over their donors
                 m=network%donors(l,k)
                 flag_not_in_list=1
                 do n=j,dep_nnode
                    if (dep_index(n).eq.m) flag_not_in_list=0
                 enddo
                 if (flag_not_in_list.eq.1.and.geometry%z(m).le.ztest(j)) then !if a donor is lower than a test elevation of its receiver, need to redistribute the sediments.
                    go to 151 !add a new node to the sedimenting node and test new elevations
                 endif
              enddo
           endif
        enddo
        go to 153 !if the loop was exited normally it means that there are no reverse gradients

151     next_node_z=maxval(geometry%z)+1.d0
        next_node_ind=-1
        index_of_receiver=-1
        do j=1,dep_nnode !add to the arrays the next node according to elevation order
           k=dep_index(j)
           if (network%ndon(k).ne.0) then
              do l=1,network%ndon(k) !loop over donors
                 m=network%donors(l,k)
                 if (geometry%z(m).lt.next_node_z)then
                    !check that it is not already in the list
                    flag_not_in_list=1
                    do n=j,dep_nnode
                       if (dep_index(n).eq.m) flag_not_in_list=0
                    enddo
                    if (flag_not_in_list.eq.1) then
                       next_node_z = geometry%z(m)
                       next_node_ind = m
                       index_of_receiver=j !index of receiver in the dep_index list
                    endif
                 endif
              enddo
           endif
        enddo
        !after the loop we have found the lowest donor. Now we can add it to the list
        !print*, 'adding node',next_node_ind,'to depositing nodes' 
        dep_nnode=dep_nnode+1
        dep_index(dep_nnode)=next_node_ind
        areas(dep_nnode)=geometry%surface(next_node_ind)
        length=dsqrt((geometry%x(next_node_ind)-geometry%x(dep_index(index_of_receiver)))**2.d0&
             &+(geometry%y(next_node_ind)-geometry%y(dep_index(index_of_receiver)))**2.d0)
        lengths(dep_nnode)=lengths(index_of_receiver) + length
        elevs(dep_nnode)=geometry%z(next_node_ind)

        sum1=depgrad*sum(areas(2:dep_nnode)*lengths(2:dep_nnode))
        sum2=sum(areas(2:dep_nnode)*elevs(2:dep_nnode))
        sum3=elevs(1)*sum(areas(2:dep_nnode))
        sum4=sum(areas(1:dep_nnode))
        dz(1)=(sedV-sum1+sum2-sum3)/sum4
!!$        if (geometry%z(dep_index(1)).eq.geometry%z(dep_index(2))) then
!!$           print*, sum1, sum2, sum3, sum4
!!$        endif
        do j=2,dep_nnode
           dz(j)=depgrad*lengths(j)-elevs(j)+elevs(1)+dz(1)
        enddo
       

        do j=1,dep_nnode
           ztest(j)=geometry%z(dep_index(j))+dz(j)
        enddo
        go to 152 !check again that the gradients have note reversed

153     sumS=0.d0 !for debuging to make sure that the correct volume is being sedimented
        do  j=1,dep_nnode
           sumS=sumS+(ztest(j)-geometry%z(dep_index(j)))*geometry%surface(dep_index(j)) 
        enddo
        !if (sumS.ne.sedV) then
        !print*, 'Available sediments volume:', sedV
        !print*, 'trying to deposit volume of:',sumS
        !endif
        !if (dep_nnode.gt.1) then
        !   print*, 'depositing', sedV, 'sediments'
        !   print*, 'over',dep_nnode,'nodes'
        !   do j=2,dep_nnode
        !      print*, 'slope of lake is',(geometry%z(dep_index(j)) - geometry%z(dep_index(1)))/lengths(j)
        !   enddo
        !endif
        do j=1,dep_nnode !here the new elevation is being updated
           geometry%z(dep_index(j)) = ztest(j)
        enddo
        deallocate(dep_index,ztest,areas,lengths,elevs,dz)

     endif
  enddo




  ! check for rivers flowing uphill after lake algorithm
  do i = 1, geometry%nnode
     if(network%receiver(i).ne.0.and.geometry%fix(i).eq.0)then
        if(geometry%z(i).lt.geometry%z(network%receiver(i)))then
           print*, 'Warning: River flowing uphill in erode.f90 after lake algorithm',i,network%receiver(i),geometry%z(i),geometry%z(network%receiver(i))
        endif
        if(geometry%z(i).eq.geometry%z(network%receiver(i)))then
           print*, 'Warning: River flowing on a plane in erode.f90 after lake algorithm'
           print*, i,network%receiver(i),geometry%z(i),geometry%z(network%receiver(i))
           print*, network%ndon(i), network%ndon(network%receiver(i)),network%receiver(network%receiver(i))
        endif
     endif
  enddo







  !
  ! Fill lakes algorithm -sdw-15.06.2010
  !    first deposit all sediment at low point of drainage basin
  !
!!$  do i=1,geometry%nnode
!!$     if(network%receiver(i).eq.0.and.geometry%fix(i).eq.0)then   !find lake nodes
!!$        sedheight=params%deltat*geometry%sediment_flux(i)/geometry%surface(i)  !maximum sed available
!!$        geometry%z(i)=geometry%z(i)+sedheight
!!$     endif
!!$  enddo

  !    then pass through stack looking for negative gradients and diffusing them until water flows down original path

!!$  lakeflag=1
!!$  icount=0
!!$  do while (lakeflag.eq.1)
!!$     icount=icount+1
!!$     lakeflag=0
!!$     do is=1,stack%nnode
!!$        i=stack%order(is)
!!$        j=network%receiver(i)
!!$        if(j.ne.0.and.geometry%fix(i).eq.0)then
!!$           if((geometry%z(j).gt.geometry%z(i)))then    !water now flowing uphill
!!$              lakeflag=1
!!$
!!$              !  correct height of both points so that sediment surface takes a specified gradient downhill
!!$
!!$              length=dsqrt((geometry%x(j)-geometry%x(i))**2+(geometry%y(j)-geometry%y(i))**2)
!!$              delz=length*depgrad/2.
!!$              zj=geometry%z(j)
!!$              zi=geometry%z(i)
!!$              areaj=geometry%surface(j)
!!$              areai=geometry%surface(i)
!!$              geometry%z(i)=(areaj*(zj+delz)+areai*(zi-delz))/(areaj+areai)+delz
!!$              geometry%z(j)=geometry%z(i)-2*delz
!!$           endif
!!$        endif
!!$     enddo
!!$     if(icount.gt.5000)then
!!$        print*, 'Warning: icount in lake fill loop exceeded'
!!$        go to 10
!!$     endif
!!$  enddo
!!$10 continue



  if (params%n.eq.1.d0) then
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).eq.0 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)
           deltal=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0)
           fact=params%deltat*geometry%k(i)/deltal*geometry%discharge(i)**params%m
           ! "erosion" is the new height of node i
           erosion=(geometry%z(i)+fact*geometry%z(j))/(1.d0+fact)
           ! the next 6 or so lines are patches to catch errors in gradient or overextrapolation of rate..
           if (geometry%z(i).lt.geometry%z(j)) then
              print*, 'error in erode.f90 -reciever higher than donor',i,j
              pause 'error in erode.f90 -reciever higher than donor'
           endif

           erosion=min(erosion,geometry%z(i))
          
           if (erosion.le.geometry%z(j)) then           
              erosion=geometry%z(j)
           endif
        
           geometry%erosion_rate(i)=-(erosion-geometry%z(i))/params%deltat

           if(geometry%erosion_rate(i).lt.0.d0)then
              geometry%erosion_rate(i)=0.d0
           else
              geometry%z(i)=erosion
           endif
        endif
     enddo
  elseif (params%n.eq.2.d0) then
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).eq.0 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)
           deltal=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0)
           fact=params%deltat*geometry%k(i)/deltal**params%n*geometry%discharge(i)**params%m
           factp=(1.d0-2.d0*fact*geometry%z(j))**2.d0-4.d0*fact*(fact*geometry%z(j)**2.d0-geometry%z(i))           
           if (factp.le.0.d0) stop 'factp negative in erode.f90'
           geometry%erosion_rate(i)=0.d0
           if (fact.ne.0.d0) then
              erosion = (2*geometry%z(j)*fact - 1.d0 + dsqrt(factp))/2.d0/fact !first root
              !erosion2 = (2*geometry%z(j)*fact - 1. - sqrt(factp))/2.d0/fact !second root  
           endif
           if (erosion .gt. geometry%z(i)) then ! I think that due to tiny numerical errors
              geometry%erosion_rate(i) = 0.0d0
              print*, 'erosion causes the node to rise'
           else
              geometry%erosion_rate(i) = (geometry%z(i) - erosion)/params%deltat
              geometry%z(i)=erosion
           endif
        endif
     enddo
  else !n > 2
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).eq.0 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)           
           deltal=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0)

           fact = params%deltat*geometry%k(i)*geometry%discharge(i)**params%m 
           erosion = rtsafe(geometry%z(i),geometry%z(j),1.0d-10,deltal,fact,params%n)
!!!! checking that all is normal
!!$           if ((params%istep/params%freq)*params%freq.eq.params%istep) then
!!$              print*, i,geometry%x(i),geometry%y(i),geometry%z(i)
!!$              print*, j,geometry%x(j),geometry%y(j),geometry%z(j)
!!$              print*, deltal,geometry%k(i),geometry%discharge(i),params%m,params%n
!!$              print*, params%deltat,erosion,(geometry%z(i) - erosion)/params%deltat
!!$              print*, '^^^^^^^^'
!!$           endif
           if (erosion.lt.geometry%z(j))then
              print*, 'source is now lower than reciever'
           endif

           if (erosion .gt. geometry%z(i)) then ! I think that due to tiny numerical errors
              geometry%erosion_rate(i) = 0.0d0
              print*, 'erosion causes the node to rise'
           else
              geometry%erosion_rate(i) = (geometry%z(i) - erosion)/params%deltat
              geometry%z(i)=erosion
           endif
        endif
     enddo
  endif

  ! set erosion rate to zero for nodes with no receiver (lakes)
  ! set the erosion rate for boundary nodes to specified value
  !   this is redundent - done above already
  do i=1, geometry%nnode
     if(geometry%fix(i).ne.0) geometry%erosion_rate(i)=min(geometry%w(i),1.d-9)
     if(geometry%erosion_rate(i).lt.1.d-9) geometry%erosion_rate(i) = 1.d-9
     if(network%receiver(i).eq.0) geometry%erosion_rate(i)=0.0d0
  enddo


 



  call time_out ('erode')

  return

end subroutine erode


!!!!!!!!!!!!!!!!!! subroutine rtsafe Newton-Raphson Method with bisection !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Based on Numerical Recipes in Fortran (page 359)     !!!!!!!!!!!!!!!!!!

function rtsafe(x1,x2,xacc,deltal,fact,n)

integer MAXIT
double precision rtsafe, x1, x2, xacc, deltal, fact,n
!EXTERNAL funcd
PARAMETER (MAXIT = 100)
integer k
double precision df,dx,dxold,f,fh,fl,temp,xh,xl


call funcd(x1,fl,df,deltal,x1,x2,fact,n)
call funcd(x2,fh,df,deltal,x1,x2,fact,n)
if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))&
stop 'root must be bracketed in rtsafe'
if (fl.eq.0.d0)then
   rtsafe = x1
   return
else if (fh.eq.0.d0)then
   rtsafe = x2
   return
else if (fl.lt.0.d0) then
   xl = x1
   xh = x2
else
   xh = x1
   xl = x2
endif
rtsafe = .5d0*(x1+x2)
dxold = abs(x2-x1)
dx=dxold
call funcd(rtsafe,f,df,deltal,x1,x2,fact,n)
do k= 1,MAXIT
   if (((rtsafe-xh)*df - f)*((rtsafe - xl)*df - f).gt.0.d0&
        .or. abs(2.d0*f).gt.abs(dxold*df)) then
      dxold = dx
      dx = 0.5d0*(xh-xl)
      rtsafe = xl+dx
      if(xl.eq.rtsafe) return
   else
      dxold=dx
      dx=f/df
      temp=rtsafe
      rtsafe = rtsafe-dx
      if (temp.eq.rtsafe) return
   endif
   if (abs(dx).lt.xacc) return
   call funcd(rtsafe,f,df,deltal,x1,x2,fact,n)
   if(f.lt.0.d0) then
      xl = rtsafe
   else
      xh = rtsafe
   endif
!   print*, 'at time',k,'guess is', rtsafe, dx, xacc
enddo
stop 'rtsafe exceeding maximum iterations'
return
end FUNCTION rtsafe


!!!!!!!!!!!! subroutine funcd !!!!!!!!!!!!!!!!!
!!! returns the erosion function and its derivative at point x

subroutine funcd (x,f,df,deltal,zi,zj,fact,n)
 
  use definitions
  implicit none
  double precision x,f,df,deltal,zi,zj,fact,n
  f = x - zi + fact*((x-zj)/deltal)**n
  df = 1.d0 + fact*n*(1/deltal)*((x-zj)/deltal)**(n-1)
  return 
end subroutine funcd
