module momentumod
contains
        ! Subrutina que calcula la energía cinética !
	real function Kinetic_E(N_atoms,Velocity_mat)
                implicit none
                !En unidades reducidas la massa es 1!
	            real, dimension(N_atoms,3), intent(in) :: Velocity_mat
                integer, intent(in) :: N_atoms
	            real :: Norma
                integer :: ii

                Kinetic_E = 0.
                do ii=1, N_atoms
                        Norma = Velocity_mat(ii,1)**2 + Velocity_mat(ii,2)**2 + Velocity_mat(ii,3)**2
                        Kinetic_E = Kinetic_E + Norma
                end do
        Kinetic_E = 0.5*Kinetic_E
        end function

        ! Subrutina para calcular el momento !
	real function P_total(N_atoms, Velocity_mat)
                implicit none
	            real, dimension(N_atoms,3), intent(in) :: Velocity_mat
                integer, intent(in) :: N_atoms
	            real :: px_total, py_total, pz_total
                integer :: ii

                px_total = 0.0
                py_total = 0.0
                pz_total = 0.0
                do ii=1, N_atoms
                        px_total = px_total + Velocity_mat(ii,1)
                        py_total = py_total + Velocity_mat(ii,2)
                        pz_total = pz_total + Velocity_mat(ii,3)
                end do
                P_total = sqrt(px_total**2 + py_total**2 + pz_total**2)
        end function
end module momentumod
