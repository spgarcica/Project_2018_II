module initial_velocitymod
use mtmod
contains
        !----------------------------------------------------------------------------------------------!
        !       Subroutine that generates normally distributed (around ± Sigma) inicial velocities     !
        !----------------------------------------------------------------------------------------------!	
        subroutine Initial_Velocity(N_atoms,num_proc,Sigma,Velocity_mat,partial_sum,Force_mat)
                implicit none
                integer, intent(in) :: N_atoms, num_proc
		real, intent(in) :: Sigma
		real, dimension(:,:), allocatable, intent(inout) :: Velocity_mat, Force_mat
		real, dimension(3), intent(out) :: partial_sum
		real, dimension(3) :: random
                integer :: ii

                allocate(Velocity_mat(N_atoms/num_proc,3))
                allocate(Force_mat(N_atoms,3))
                Velocity_mat = 1
                do ii=1, N_atoms/num_proc
                        random(1) = real(grnd())
                        random(2) = real(grnd())
                        random(3) = real(grnd())
                        Velocity_mat(ii,:) = Velocity_mat(ii,:) + random*Sigma
                end do
                do ii=1, 3
                        partial_sum(ii) = sum(Velocity_mat(:,ii))
                end do
        end subroutine
end module initial_velocitymod
