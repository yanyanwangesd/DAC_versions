subroutine initialize_parameters (params)

! This routine is used to initialize the model parameters
! i.e. those that do not have a spatial variation (I;E. not attached to a node)
! All of these parameters are stored in a (parm) derived type params

use definitions

implicit none

type (parm) params

double precision u,p,slope,k,sloped,area, sloped_total
integer nsteps, nout

call time_in ('initialize_parameters')
   

!  echo parameters to screen
!

print*, ''
print*, 'DAC (Divide And Capture). Code branch with polynomial input for different fields'
print*, ''
print*, ''
print*, 'Timestepping parameters:'
print*, 'delta t =', params%deltat
print*, 'endtime = ',params%tfinal
print*,  'frequency of output =', params%freq
nsteps=params%tfinal/params%deltat
nout=nsteps/params%freq+1
print*, nsteps, ' timesteps with ',nout,' outputs'

print*, ''
print*, 'Fluvial erosion parameters: '
print*, 'n = ',params%n
print*, 'm = ',params%m
print*, 'Hacks law h = ',params%h
print*, 'Hacks Law ka = ',params%ka
print*, 'K = ',params%k_scalar1
print*, 'Precip = ',params%rainfall_height

print*, ''
print*, 'Channel head erosion parameters: '
print*, 'Xc = ',params%xc
print*, 'tan theta = ',params%tanthetac
print*, 'Diffusivity = ',params%diffusivity/2.


print*, ''
print*, 'Uplift Parameters: '
print*, 'Interior Uplift rate = ',params%uplift_scalar1
print*, 'Boundary Uplift rate = ',params%uplift_scalar2
print*, ''
print*, 'Adding and Removing Node parameters: '
print*, 'Maximum allowable channel length = ',params%lmax
print*, 'Minimum allowable channel head area = ',params%amin
print*, 'Max Divide length = ',params%ldivmax
print*, ''
print*, ''
print*, ''
! calculate at what length the channel slope will become steeper than the imposed channel head
! this is a function of incision rate which is specified in subroutine uplift_advect and precip given in find_precip
! and k given in initialize_geometry. 
! Give them again here as a temp variables just for this calc. This is done in SS here.

u=params%uplift_scalar1
p=params%rainfall_height
k=params%k_scalar1
slope=(u/(k*params%xc**(params%h*params%m)*p**params%m*params%ka**params%m))**(1./params%n)
slope=180*atan(slope)/3.141592
sloped=u*params%xc/(params%diffusivity/2.)
sloped_total=u*params%xc/(params%diffusivity)
sloped=180*atan(sloped)/3.141592
sloped_total=180*atan(sloped_total)/3.141592
if(slope.gt.(180*atan(params%tanthetac)/3.141592))then
print*, '!! warning channel head slope is less steep than channel'
endif
if(sloped_total.lt.slope)then
print*, 'Warning!!! diffusive channel head slope is less steep than channel'
endif
print*, 'max slope of channel: ', slope
print*, 'slope of channel head: ', (180*atan(params%tanthetac)/3.141592)
print*, 'max diffusive slope of channel head: ', sloped_total

! check the channel head area

area=params%ka*params%xc**params%h
if(params%amin.lt.area)print*,'Warning!!!! minimum channel head area too small'
print*, 'Channel Head Area = ',area/1.d6
print*, 'Remesh minimum area used = ',params%amin/1.d6

!check linearity when using transient divide solution
if (params%transient_divide.and.params%n.ne.1) then
   print*, 'Cant solve transient divide with n=',params%n
   print*, 'Use slope exponent n=1'
endif


call time_out ('initialize_parameters')

return

end subroutine initialize_parameters
