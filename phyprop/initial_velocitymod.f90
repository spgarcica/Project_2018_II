module initial_velocitymod
use mtmod
contains
        !----------------------------------------------------------------------------------------------!
        !       Subroutine that generates normally distributed (around ± Sigma) inicial velocities     !
        !----------------------------------------------------------------------------------------------!	
        subroutine Initial_Velocity(N_atoms,Sigma,Velocity_mat)
                implicit none
                integer, intent(in) :: N_atoms
                real, intent(in) :: Sigma
                real, dimension(N_atoms,3), intent(out) :: Velocity_mat
                real, dimension(3) :: random
                integer :: ii

               Velocity_mat = 1
               do ii=1, N_atoms
                       random(1) = real(grnd())
                       random(2) = real(grnd())
                       random(3) = real(grnd())
                       Velocity_mat(ii,:) = Velocity_mat(ii,:) + random * Sigma
               end do
               do ii=1, 3
                       Velocity_mat(:,ii) = Velocity_mat(:,ii) - sum(Velocity_mat(:,ii))/N_atoms
               end do
        end subroutine
end module initial_velocitymod
