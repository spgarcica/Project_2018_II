module gaussmod
contains
        ! Realizado -> Alejandro !
        ! Función de la gaussiana !
        real function gaussiana (Temp, x)
                implicit none
                real, intent(in) :: Temp, x
                real :: pi
                parameter (pi = 4.D0*DATAN(1.D0))

                gaussiana = (1.0/sqrt(2.0*pi*Temp))*exp(-(x**2)/(2.0*Temp))
        end function gaussiana
        
        ! Realizado -> Alejandro !
        ! Sencilla subrutina para crear un archivo con una gaussiana !
        ! de manera analítica y comparar con las velocidades         !
        subroutine gauss (Temp, N_atoms, NSteps)
                implicit none
                real, intent(in) :: Temp
                integer, intent(in) :: N_atoms, NSteps
                real :: x, incr_x
                integer :: i, N

                incr_x = 0.01
                N=20.0/incr_x

                open(unit=25, file='gauss.dat',form='formatted',status='unknown')
                write(25,*) "#Formato: x, f(x)"
                do i=0, N
                        x = -10 + i*incr_x
                        write(25,*) x, (N_atoms*NSteps*gaussiana(Temp,x))/5.0
                end do
                close(25)
        end subroutine
end module gaussmod
