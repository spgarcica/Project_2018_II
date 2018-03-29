module AndersenMod
   use normaldist 
   use mtmod
   implicit none
   contains
        !----------------------------------------------------------------------------!
        !                        Andersen thermostat subroutine                      !
        !----------------------------------------------------------------------------!
        subroutine Andersen(N_atoms,myrank,num_proc,A_Prob,Vel_mat,Tinst)
                implicit none
                real, dimension(N_atoms/num_proc,3), intent(inout) :: Vel_mat
                integer, intent(in) :: N_atoms, num_proc, myrank
                real, intent(in) :: Tinst, A_Prob
                real, dimension(3) :: probvec
                integer :: ii, jj
                integer :: partition
                real :: random, sigma23
                partition = N_atoms/num_proc
                jj = 0

                if (A_prob > 0) then
                ! Imposing the system to have a Tins !
                sigma23 = sqrt(Tinst)
                do ii= (partition*myrank)+1,(myrank+1)*partition
                        jj = jj + 1
                        Random = real(grnd())
                        if (random .lt. A_Prob) then
                                call r4vec_normal_ab(3,0.,Sigma23,probvec)
                                Vel_mat(jj,:) = probvec
                        end if
                end do
                end if
        end subroutine

end module
