module momentummod
contains
        real function Kinetic_E(N_atoms,Vel_mat_total,rank,param)
        implicit none
        integer, intent(in) :: N_atoms, rank, param
        real :: Norma
        integer :: ii
        real, dimension(N_atoms,3), intent(in) :: Vel_mat_total

        Kinetic_E = 0.0
        do ii=(rank*param)+1, (rank+1)*param
                Norma = Vel_mat_total(ii,1)**2 + Vel_mat_total(ii,2)**2 + Vel_mat_total(ii,3)**2
                Kinetic_E = Kinetic_E + Norma
        end do
        Kinetic_E = 0.5*Kinetic_E

        end function

        subroutine P_sum(N_atoms,Vel_mat_total,rank,param,partial_P)
        implicit none
        real, dimension(N_atoms,3), intent(in) :: Vel_mat_total
        real, dimension(3), intent(inout) :: partial_P
        integer, intent(in) :: N_atoms, rank, param
        integer :: ii

        partial_P = 0.0
        do ii=(rank*param)+1, (rank+1)*param
                partial_P(1) = partial_P(1) + Vel_mat_total(ii,1)
                partial_P(2) = partial_P(2) + Vel_mat_total(ii,2)
                partial_P(3) = partial_P(3) + Vel_mat_total(ii,3)
        end do
        end subroutine
end module momentummod
