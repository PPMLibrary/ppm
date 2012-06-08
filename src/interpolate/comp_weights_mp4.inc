!Compute kernel weights


!-------------------------
! Long ugly version
! In case the mesh array does not have a halo layer, e.g. with patches, we
! cannot apply a centered remeshing stencil everywhere. Rather than treating the
! borders separately, we hack the code so that out-of-bound indices
! are shifted back into the allocated region while their kernel weights are set
! to zero. This way, instead of accessing node (-1,1,1), we access node (1,1,1)
! but set a10 to zero. Since the weight is a10a20a30=a10*a20*a30=0, the
! contribution of this node will be nil.
! TODO (try replacing the homemade Heaviside functions
! by MIN or MAX statements. That probably also vectorizes and may
! well be faster...)
! Alternatively, one can use the shorter version (below) and discard the
! particles that are near the patch boundaries. This should be done by 
! replacing all occurences of p%end_ext in part/part_interp_to_mesh.f by
! p%end_red   (instead of selecting particles that are in the extended patch
! domain, we select these that are in the reduced patch, ie patch minus
! ghostsize)
!-------------------------

!--------------
!  x component
!--------------
x01 = xp(1,iq)*dxxi-p%istart(1) + 1

ip10 = FLOOR(x01)
ip11 = ip10 + 1
ip12 = ip11 + 1
ip13 = ip11 + 2

ip10l = ip10 - p%lo_a(1)
ip10h = p%hi_a(1) - ip10

! some kind of Heaviside function
!    for particles near the left boundary

l10 = (1 + SIGN(one,ip10l  )) / 2
l11 = (1 + SIGN(one,ip10l+1)) / 2
l12 = (1 + SIGN(one,ip10l+2)) / 2

!    for particles near the right boundary
h11 = (1 + SIGN(one,ip10h-1)) / 2
h12 = (1 + SIGN(one,ip10h-2)) / 2
h13 = (1 + SIGN(one,ip10h-3)) / 2

x10 = (x01-REAL(ip10,mk)) * l10       + 1._mk
x11 = (x01-REAL(ip11,mk)) * l11 * h11 + 1._mk
x12 = (x01-REAL(ip10,mk)) * l12 * h12 - 1._mk
x13 = (x01-REAL(ip11,mk))       * h13 - 1._mk

ip10 = (ip10-p%lo_a(1) + ABS(ip10-p%lo_a(1)))/2 + p%lo_a(1)
ip11 = (ip11-p%lo_a(1) + ABS(ip11-p%lo_a(1)))/2 + p%lo_a(1)
ip12 = (ip12-p%lo_a(1) + ABS(ip12-p%lo_a(1)))/2 + p%lo_a(1)

ip11 = p%hi_a(1) - (p%hi_a(1)-ip11 + ABS(p%hi_a(1)-ip11))/2
ip12 = p%hi_a(1) - (p%hi_a(1)-ip12 + ABS(p%hi_a(1)-ip12))/2
ip13 = p%hi_a(1) - (p%hi_a(1)-ip13 + ABS(p%hi_a(1)-ip13))/2

!--------------
!  y component
!--------------
x02 = xp(2,iq)*dxyi-p%istart(2) + 1

ip20 = FLOOR(x02)
ip21 = ip20 + 1
ip22 = ip21 + 1
ip23 = ip21 + 2

ip20l = ip20 - p%lo_a(2)
ip20h = p%hi_a(2) - ip20

! some kind of Heaviside function
!    for particles near the left boundary
l20 = (1 + SIGN(one,ip20l  )) / 2
l21 = (1 + SIGN(one,ip20l+1)) / 2
l22 = (1 + SIGN(one,ip20l+2)) / 2

!    for particles near the right boundary
h21 = (1 + SIGN(one,ip20h-1)) / 2
h22 = (1 + SIGN(one,ip20h-2)) / 2
h23 = (1 + SIGN(one,ip20h-3)) / 2

x20 = (x02-REAL(ip20,mk)) * l20       + 1._mk
x21 = (x02-REAL(ip21,mk)) * l21 * h21 + 1._mk
x22 = (x02-REAL(ip20,mk)) * l22 * h22 - 1._mk
x23 = (x02-REAL(ip21,mk))       * h23 - 1._mk

ip20 = (ip20-p%lo_a(2) + ABS(ip20-p%lo_a(2)))/2 + p%lo_a(2)
ip21 = (ip21-p%lo_a(2) + ABS(ip21-p%lo_a(2)))/2 + p%lo_a(2)
ip22 = (ip22-p%lo_a(2) + ABS(ip22-p%lo_a(2)))/2 + p%lo_a(2)

ip21 = p%hi_a(2) - (p%hi_a(2)-ip21 + ABS(p%hi_a(2)-ip21))/2
ip22 = p%hi_a(2) - (p%hi_a(2)-ip22 + ABS(p%hi_a(2)-ip22))/2
ip23 = p%hi_a(2) - (p%hi_a(2)-ip23 + ABS(p%hi_a(2)-ip23))/2


!--------------
!  z component
!--------------

x03 = xp(3,iq)*dxzi-p%istart(3) + 1

ip30 = FLOOR(x03)
ip31 = ip30 + 1
ip32 = ip31 + 1
ip33 = ip31 + 2

ip03l = ip30 - p%lo_a(3)
ip03h = p%hi_a(3) - ip30

! some kind of Heaviside function
!    for particles near the left boundary
l30 = (1 + SIGN(one,ip03l  )) / 2
l31 = (1 + SIGN(one,ip03l+1)) / 2
l32 = (1 + SIGN(one,ip03l+2)) / 2

!    for particles near the right boundary
h31 = (1 + SIGN(one,ip03h-1)) / 2
h32 = (1 + SIGN(one,ip03h-2)) / 2
h33 = (1 + SIGN(one,ip03h-3)) / 2

x30 = (x03-REAL(ip30,mk)) * l30       + 1._mk
x31 = (x03-REAL(ip31,mk)) * l31 * h31 + 1._mk
x32 = (x03-REAL(ip30,mk)) * l32 * h32 - 1._mk
x33 = (x03-REAL(ip31,mk))       * h33 - 1._mk

ip30 = (ip30-p%lo_a(3) + ABS(ip30-p%lo_a(3)))/2 + p%lo_a(3)
ip31 = (ip31-p%lo_a(3) + ABS(ip31-p%lo_a(3)))/2 + p%lo_a(3)
ip32 = (ip32-p%lo_a(3) + ABS(ip32-p%lo_a(3)))/2 + p%lo_a(3)

ip31 = p%hi_a(3) - (p%hi_a(3)-ip31 + ABS(p%hi_a(3)-ip31))/2
ip32 = p%hi_a(3) - (p%hi_a(3)-ip32 + ABS(p%hi_a(3)-ip32))/2
ip33 = p%hi_a(3) - (p%hi_a(3)-ip33 + ABS(p%hi_a(3)-ip33))/2

!-------------------------
! Shorter version
! (but leads to out-of-bound access if the array is not padded with a
! halo layer. This is the case for a patch that does not cover the whole
! domain, for example).
!-------------------------
!                 x01 = xp(1,iq)*dxxi-p%istart(1) + 1
!                 x02 = xp(2,iq)*dxyi-p%istart(2) + 1
!                 x03 = xp(3,iq)*dxzi-p%istart(3) + 1
!
!                 ip10 = FLOOR(x01)
!                 ip20 = FLOOR(x02)
!                 ip30 = FLOOR(x03)
!
!
!                 ip11 = ip10 + 1
!                 ip21 = ip20 + 1
!                 ip31 = ip30 + 1
!
!                 ip12 = ip11 + 1
!                 ip22 = ip21 + 1
!                 ip32 = ip31 + 1
!
!                 ip13 = ip11 + 2
!                 ip23 = ip21 + 2
!                 ip33 = ip31 + 2
!
!                 xp1 = x01-REAL(ip10,mk)
!                 xp2 = x02-REAL(ip20,mk)
!                 xp3 = x03-REAL(ip30,mk)
!
!                 x10 = xp1 + 1.0_mk
!                 x11 = x10 - 1.0_mk
!                 x12 = x10 - 2.0_mk
!                 x13 = x10 - 3.0_mk
!
!                 x20 = xp2 + 1.0_mk
!                 x21 = x20 - 1.0_mk
!                 x22 = x20 - 2.0_mk
!                 x23 = x20 - 3.0_mk
!
!                 x30 = xp3 + 1.0_mk
!                 x31 = x30 - 1.0_mk
!                 x32 = x30 - 2.0_mk
!                 x33 = x30 - 3.0_mk

!----------
! weights of the M4p kernel
!----------

a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

a10a20a30 = a10*a20*a30
a10a20a31 = a10*a20*a31
a10a20a32 = a10*a20*a32
a10a20a33 = a10*a20*a33
a10a21a30 = a10*a21*a30
a10a21a31 = a10*a21*a31
a10a21a32 = a10*a21*a32
a10a21a33 = a10*a21*a33
a10a22a30 = a10*a22*a30
a10a22a31 = a10*a22*a31
a10a22a32 = a10*a22*a32
a10a22a33 = a10*a22*a33
a10a23a30 = a10*a23*a30
a10a23a31 = a10*a23*a31
a10a23a32 = a10*a23*a32
a10a23a33 = a10*a23*a33
a11a20a30 = a11*a20*a30
a11a20a31 = a11*a20*a31
a11a20a32 = a11*a20*a32
a11a20a33 = a11*a20*a33
a11a21a30 = a11*a21*a30
a11a21a31 = a11*a21*a31
a11a21a32 = a11*a21*a32
a11a21a33 = a11*a21*a33
a11a22a30 = a11*a22*a30
a11a22a31 = a11*a22*a31
a11a22a32 = a11*a22*a32
a11a22a33 = a11*a22*a33
a11a23a30 = a11*a23*a30
a11a23a31 = a11*a23*a31
a11a23a32 = a11*a23*a32
a11a23a33 = a11*a23*a33
a12a20a30 = a12*a20*a30
a12a20a31 = a12*a20*a31
a12a20a32 = a12*a20*a32
a12a20a33 = a12*a20*a33
a12a21a30 = a12*a21*a30
a12a21a31 = a12*a21*a31
a12a21a32 = a12*a21*a32
a12a21a33 = a12*a21*a33
a12a22a30 = a12*a22*a30
a12a22a31 = a12*a22*a31
a12a22a32 = a12*a22*a32
a12a22a33 = a12*a22*a33
a12a23a30 = a12*a23*a30
a12a23a31 = a12*a23*a31
a12a23a32 = a12*a23*a32
a12a23a33 = a12*a23*a33
a13a20a30 = a13*a20*a30
a13a20a31 = a13*a20*a31
a13a20a32 = a13*a20*a32
a13a20a33 = a13*a20*a33
a13a21a30 = a13*a21*a30
a13a21a31 = a13*a21*a31
a13a21a32 = a13*a21*a32
a13a21a33 = a13*a21*a33
a13a22a30 = a13*a22*a30
a13a22a31 = a13*a22*a31
a13a22a32 = a13*a22*a32
a13a22a33 = a13*a22*a33
a13a23a30 = a13*a23*a30
a13a23a31 = a13*a23*a31
a13a23a32 = a13*a23*a32
a13a23a33 = a13*a23*a33

