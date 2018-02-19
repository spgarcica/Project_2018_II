MODULE pbc
   implicit none

   contains
        !--------------------------------------------------------------------------!
        ! Function to correct the pbc subroutine when the length the atom has      !
        ! moved is higher than L.                                                  !
        !--------------------------------------------------------------------------!
        real function PBC_Cor(Val,Siz)
                implicit none
                real, intent(in) :: Val, Siz
                PBC_Cor = abs(aint(Val/Siz)) * Siz
        end function PBC_Cor

        !--------------------------------------------------------------------------!
        !         Subroutine that applies the periodic boundary conditions         !
        !--------------------------------------------------------------------------!
        subroutine pbc(coord, N, L) !No sé si debo poner el out coord2
                implicit none
                real, dimension(N,3), intent(inout) :: coord
                real, intent(in) :: L
                integer, intent(in) :: N
                integer :: jN
                real, dimension(N) :: x, y, z

                !Separa la matriz en 3 vectores
                x = coord(:,1)
                y = coord(:,2)
                z = coord(:,3)

                do jN = 1, N
                        if (x(jN) .lt. 0) then
                                x(jN) = x(jN) + PBC_Cor(x(jN),L) + L
                        else if (x(jN) .gt. L) then
                                x(jN) = x(jN) - PBC_Cor(x(jN),L)
                        end if
                        
                        if (y(jN) .lt. 0) then
                                y(jN) = y(jN) + PBC_Cor(y(jN),L) + L
                        else if (y(jN) .gt. L) then
                                y(jN) = y(jN) - PBC_Cor(y(jN),L)
                        end if

                        if (z(jN) .lt. 0) then
                                z(jN) = z(jN) + PBC_Cor(z(jN),L) + L
                        else if (z(jN) .gt. L) then
                               z(jN) = z(jN) - PBC_Cor(z(jN),L)
                        end if
                end do

                !Ponemos las condiciones 
                coord(:,1) = x
                coord(:,2) = y
                coord(:,3) = z
        end subroutine
END MODULE pbc
