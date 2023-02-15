! translated from http://www.skycoyote.com/cntr/


subroutine calculate_centroid (xy, ctr)

real*8 xy(2,3), ctr(2)
real*8 x1, x2, x3, x4, y1, y2, y3, y4
real*8 a1, a2, b1, b2

x1 = xy(1,1)	! get endpoints of first median
y1 = xy(2,1)	! get endpoints of first median
x2 = (xy(1,2) + xy(1,3))/2. ! get endpoints of first median
y2 = (xy(2,2) + xy(2,3))/2. ! get endpoints of first median

x3 = xy(1,2)	! get endpoints of second median  (only need two)
y3 = xy(2,2)	! get endpoints of second median
x4 = (xy(1,3) + xy(1,1))/2. ! get endpoints of second median
y4 = (xy(2,3) + xy(2,1))/2. ! get endpoints of second median

! see if either median is vertical (slope == infinity)
if (x1.eq.x2) then
   x1 = xy(1,3)	! use third median
   y1 = xy(2,3)	!  use third median
   x2 = (xy(1,2) + xy(1,1))/2. !  use third median
   y2 = (xy(2,2) + xy(2,1))/2. !  use third median
elseif (x3.eq.x4) then
   x3 = xy(1,3)	! use third median
   y3 = xy(2,3)	!  use third median
   x4 = (xy(1,2) + xy(1,1))/2. !  use third median
   y4 = (xy(2,2) + xy(2,1))/2. !  use third median
endif

a1 = (y2 - y1) / (x2 - x1)	! compute slope of first median
b1 = y1 - a1 * x1	! compute intercept of first median

a2 = (y4 - y3) / (x4 - x3)	! compute slope of second median
b2 = y3 - a2 * x3	! compute intercept of second median

! solve a1 * x + b1 = a2 * x + b2

ctr(1) = (b2 - b1) / (a1 - a2)	! solve for x coordinate of intersection
ctr(2) = a1 * ctr(1) + b1	! solve for y coordinate of intersection

return

end subroutine calculate_centroid

