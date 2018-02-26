module normaldist
        use mtmod
        contains
subroutine r4vec_normal_ab ( n, a, b, x )

!*****************************************************************************80
!
!! R4VEC_NORMAL_AB returns a scaled pseudonormal R4VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    An R4VEC is a vector of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input, real ( kind = 4 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 4 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) m, counti
  real ( kind = 4 ) r(n+1)
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = real(grnd())

    if ( r(1) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4VEC_NORMAL_AB - Fatal error!'
      write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = real(grnd())

    x(x_hi_index) = &
      sqrt ( - 2.0E+00 * log ( r(1) ) ) * cos ( 2.0E+00 * r4_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    do counti=1, 2*m
    r(counti) = real(grnd())
    end do

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1


    do counti=1, 2*m
    r(counti) = real(grnd())
    end do

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0E+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
end module
