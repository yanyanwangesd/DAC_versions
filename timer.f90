subroutine time_in (name)

! Utility routines to compute the time spent in each subroutine
!
! For this to work, each routine must start with a call to time_in and end with a call to time_out.
! The only argument to these two routines is the name of the routine itself (as a string variable).
! For example for the routine called "erode", one must call:
!      subroutine erode (geometry, network,stack,params)
!      ...
!      call time_in ('erode')
!      ...
!      call time_out ('erode')
!      return
!      end subroutine erode
!
! At the start of the main program, the first executabme line must be acall to start_timing
! and at the end of the main program one must add a call to show_timing:
!      program DivideAndCapture
!      ...
!      call start_timing
!      ...
!      call show_timing
!      end program DivideAndCapture

use definitions

implicit none

type (timr) timer
character*(*) name
integer len,i
real time0

common /global_timer/ timer

len=len_trim (name)
call cpu_time (time0)

timer%timein=time0
timer%namein(1:len)=name(1:len)

!print*,name(1:len)

  do i=1,timer%ntimer
  if (name(1:len).eq.timer%name(i)(1:len)) return
  enddo

timer%ntimer=timer%ntimer+1
if (timer%ntimer.gt.1024) stop 'too many timers...'
timer%name(timer%ntimer)(1:len)=name(1:len)

return

end subroutine time_in

!---

subroutine time_out (name)

use definitions

implicit none

type (timr) timer
character*(*) name
integer len,i
real time0

common /global_timer/ timer

len=len_trim (name)

call cpu_time (time0)

  do i=1,timer%ntimer
    if (name(1:len).eq.timer%name(i)(1:len)) then
    timer%time_spent(i)=timer%time_spent(i)+time0-timer%timein
    return
    endif
  enddo

print*,'Cannot find routine name in time_out called by '//name(1:len)
stop

end subroutine time_out

!---

subroutine start_timing

use definitions

implicit none

type (timr) timer

call cpu_time (timer%time_start)

end subroutine start_timing

!---

subroutine show_timing

use definitions

implicit none

type (timr) timer
integer i

common /global_timer/ timer

print*,'-----------------------------------------------'
print*,'TIMER INFORMATION'
print*,'-----------------------------------------------'
  do i=1,timer%ntimer
  write (*,'(a25,f20.10,a3)') timer%name(i)(1:len_trim(timer%name(i))),timer%time_spent(i),'sec'
  enddo
print*,'-----------------------------------------------'
call cpu_time (timer%time_stop)
print*,'Total time:',timer%time_stop-timer%time_start

return

end subroutine show_timing
