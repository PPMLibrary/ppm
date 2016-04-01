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
               !
               ! TODO (try replacing the homemade Heaviside functions
               ! by MIN or MAX statements. That probably also vectorizes and may
               ! well be faster...)
               ! Yaser : I have replaced the Heaviside function by a simple IF
               ! which I think is much faster.
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
               ! for particles near the left boundary

               !l10 = (1 + SIGN(1,ip10l  )) / 2
               !l11 = (1 + SIGN(1,ip10l+1)) / 2
               !l12 = (1 + SIGN(1,ip10l+2)) / 2

               IF (ip10l.LT.-2) THEN
                  l10=0
                  l11=0
                  l12=0
               ELSE IF (ip10l.LT.-1) THEN
                  l10=0
                  l11=0
                  l12=1
               ELSE IF (ip10l.LT.0) THEN
                  l10=0
                  l11=1
                  l12=1
               ELSE
                  l10=1
                  l11=1
                  l12=1
               ENDIF

               ! for particles near the right boundary
               !h11 = (1 + SIGN(1,ip10h-1)) / 2
               !h12 = (1 + SIGN(1,ip10h-2)) / 2
               !h13 = (1 + SIGN(1,ip10h-3)) / 2

               IF (ip10h.LT.1) THEN
                  h11=0
                  h12=0
                  h13=0
               ELSE IF (ip10h.LT.2) THEN
                  h11=1
                  h12=0
                  h13=0
               ELSE IF (ip10h.LT.3) THEN
                  h11=1
                  h12=1
                  h13=0
               ELSE
                  h11=1
                  h12=1
                  h13=1
               ENDIF

               x10 = (x01-REAL(ip10,MK)) * l10       + one
               x11 = (x01-REAL(ip11,MK)) * l11 * h11 + one
               x12 = (x01-REAL(ip10,MK)) * l12 * h12 - one
               x13 = (x01-REAL(ip11,MK))       * h13 - one

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
               ! for particles near the left boundary
               !l20 = (1 + SIGN(1,ip20l  )) / 2
               !l21 = (1 + SIGN(1,ip20l+1)) / 2
               !l22 = (1 + SIGN(1,ip20l+2)) / 2

               IF (ip20l.LT.-2) THEN
                  l20=0
                  l21=0
                  l22=0
               ELSE IF (ip20l.LT.-1) THEN
                  l20=0
                  l21=0
                  l22=1
               ELSE IF (ip20l.LT.0) THEN
                  l20=0
                  l21=1
                  l22=1
               ELSE
                  l20=1
                  l21=1
                  l22=1
               ENDIF

               !    for particles near the right boundary
               !h21 = (1 + SIGN(1,ip20h-1)) / 2
               !h22 = (1 + SIGN(1,ip20h-2)) / 2
               !h23 = (1 + SIGN(1,ip20h-3)) / 2

               IF (ip20h.LT.1) THEN
                  h21=0
                  h22=0
                  h23=0
               ELSE IF (ip20h.LT.2) THEN
                  h21=1
                  h22=0
                  h23=0
               ELSE IF (ip20h.LT.3) THEN
                  h21=1
                  h22=1
                  h23=0
               ELSE
                  h21=1
                  h22=1
                  h23=1
               ENDIF


               x20 = (x02-REAL(ip20,MK)) * l20       + one
               x21 = (x02-REAL(ip21,MK)) * l21 * h21 + one
               x22 = (x02-REAL(ip20,MK)) * l22 * h22 - one
               x23 = (x02-REAL(ip21,MK))       * h23 - one

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
               !l30 = (1 + SIGN(1,ip03l  )) / 2
               !l31 = (1 + SIGN(1,ip03l+1)) / 2
               !l32 = (1 + SIGN(1,ip03l+2)) / 2

               IF (ip03l.LT.-2) THEN
                  l30=0
                  l31=0
                  l32=0
               ELSE IF (ip03l.LT.-1) THEN
                  l30=0
                  l31=0
                  l32=1
               ELSE IF (ip03l.LT.0) THEN
                  l30=0
                  l31=1
                  l32=1
               ELSE
                  l30=1
                  l31=1
                  l32=1
               ENDIF

               ! for particles near the right boundary
               !h31 = (1 + SIGN(1,ip03h-1)) / 2
               !h32 = (1 + SIGN(1,ip03h-2)) / 2
               !h33 = (1 + SIGN(1,ip03h-3)) / 2

               IF (ip03h.LT.1) THEN
                  h31=0
                  h32=0
                  h33=0
               ELSE IF (ip03h.LT.2) THEN
                  h31=1
                  h32=0
                  h33=0
               ELSE IF (ip03h.LT.3) THEN
                  h31=1
                  h32=1
                  h33=0
               ELSE
                  h31=1
                  h32=1
                  h33=1
               ENDIF

               x30 = (x03-REAL(ip30,MK)) * l30       + one
               x31 = (x03-REAL(ip31,MK)) * l31 * h31 + one
               x32 = (x03-REAL(ip30,MK)) * l32 * h32 - one
               x33 = (x03-REAL(ip31,MK))       * h33 - one

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
               !                 xp1 = x01-REAL(ip10,MK)
               !                 xp2 = x02-REAL(ip20,MK)
               !                 xp3 = x03-REAL(ip30,MK)
               !
               !                 x10 = xp1 + one
               !                 x11 = x10 - one
               !                 x12 = x10 - two
               !                 x13 = x10 - three
               !
               !                 x20 = xp2 + one
               !                 x21 = x20 - one
               !                 x22 = x20 - two
               !                 x23 = x20 - three
               !
               !                 x30 = xp3 + one
               !                 x31 = x30 - one
               !                 x32 = x30 - two
               !                 x33 = x30 - three

               !----------
               ! weights of the M4p kernel
               !----------

               a10 = two + (-four+(two5-half*x10)*x10)*x10
               a20 = two + (-four+(two5-half*x20)*x20)*x20
               a30 = two + (-four+(two5-half*x30)*x30)*x30

               a11 = one + (-two5+one5*x11)*x11**2
               a21 = one + (-two5+one5*x21)*x21**2
               a31 = one + (-two5+one5*x31)*x31**2

               a12 = one + (-two5-one5*x12)*x12**2
               a22 = one + (-two5-one5*x22)*x22**2
               a32 = one + (-two5-one5*x32)*x32**2

               a13 = two + (four + (two5+half*x13)*x13)*x13
               a23 = two + (four + (two5+half*x23)*x23)*x23
               a33 = two + (four + (two5+half*x33)*x33)*x33

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

