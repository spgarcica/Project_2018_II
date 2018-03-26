module velocity_verlet2mod
use pbcmod
contains
        subroutine Velocity_Verlet2(N_atoms,myrank,num_proc,dt,Vel_mat_total,Vel_mat,Force_mat_total,NewForce_mat_total)
        implicit none
                integer :: ii, jj
                integer :: partition
                integer, intent(in) :: N_atoms, myrank, num_proc
                real, intent(IN) :: dt
               
                real, dimension(N_atoms,3), intent(in) :: Vel_mat_total, Force_mat_total, NewForce_mat_total
                real, dimension(N_atoms/num_proc,3), intent(inout) :: Vel_mat
                partition = N_atoms/num_proc
                jj = 0
                do ii=(partition*myrank)+1, (myrank+1)*partition
                         jj = jj + 1
                         Vel_mat(jj,:) = Vel_mat_total(ii,:) + 0.5*(Force_mat_total(ii,:) + NewForce_mat_total(ii,:))*dt
                end do
        end subroutine
end module velocity_verlet2mod
                
