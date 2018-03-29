module eulermod
use pbcmod
   contains
        !-----------------------------------------------------------------------------------------!
        !                         Implementation of the Euler algorithm                           !
        !-----------------------------------------------------------------------------------------!
        subroutine euler(N_atoms,myrank,num_proc,dt,Pos_mat_total,Pos_mat,Vel_mat_total,Vel_mat,Force_mat_total,L,Partition)
        implicit none
                integer :: ii, jj
                integer, intent(in) :: N_atoms, myrank, num_proc
                integer, intent(inout) :: Partition
                real, intent(in) :: dt, L
                real, dimension(N_atoms,3),intent(in) :: Pos_mat_total,Vel_mat_total,Force_mat_total
                real, dimension(N_atoms/num_proc,3),intent(inout) :: Pos_mat, Vel_mat

                jj = 0
                  
                do ii = (partition*myrank)+1, (myrank+1)*partition
                          jj = jj + 1
                          Pos_mat(jj,:) = Pos_mat_total(ii,:) + Vel_mat_total(ii,:)*dt + 0.5*Force_mat_total(ii,:)*dt**2
                          Vel_mat(jj,:) = Vel_mat_total(ii,:) + Force_mat_total(ii,:)*dt
                end do
                call pbc(Pos_mat,partition,L)
         end subroutine
end module eulermod
