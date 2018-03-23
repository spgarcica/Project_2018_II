module velocity_verletmod
use pbcmod
contains
        subroutine Velocity_Verlet(N_atoms,myrank,num_proc,dt,Pos_mat_total,Pos_mat,Vel_mat_total,Vel_mat,Force_mat_total,Force_mat,L)
        implicit none
                integer :: ii, jj
                integer :: partition
                integer, intent(in) :: N_atoms, myrank, num_proc
                real, intent(in) :: dt, L

                real, dimension(N_atoms,3), intent(in) :: Pos_mat_total, Vel_mat_total, Force_mat_total
                real, dimension(N_atoms/num_proc,3), intent(inout) :: Pos_mat, Vel_mat, Force_mat
                partition = N_atoms/num_proc
                jj=0
                do ii=(partition*myrank)+1, (myrank+1)*partition
                        jj = jj + 1
                        Pos_mat(jj,:) = Pos_mat_total(ii,:) + Vel_mat_total(ii,:)*dt + 0.5*Force_mat_total(ii,:)*dt**2
                end do
                call pbc(Pos_mat,partition,L)
        end subroutine
end module velocity_verletmod
