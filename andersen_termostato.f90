module andersen
   implicit none

   contains
        ! Realizado -> Víctor, Revisado -> Sergio !
        ! Termostato de Andersen !
        subroutine Andersen(N_atoms,A_Prob,Velocities_mat,Tinst)
                implicit none
                real, dimension(N_atoms,3), intent(inout) :: Velocities_mat
                integer, intent(in) :: N_atoms
                real, intent(in) :: Tinst, A_Prob
                real, dimension(3) :: probvec
                integer :: ii
                real :: random, sigma23

                if (A_prob > 0) then
                !Forzamos que nuestro sistema tenga Tinst.
                sigma23 = sqrt(Tinst)
                do ii=1, N_atoms !Hacemos un loop sobre todas las partículas
                        Random = real(grnd())
                        if (random .lt. A_Prob) then
                                call r4vec_normal_ab(3,0.,Sigma23,probvec)
                                Velocities_mat(ii,:) = probvec
                        end if
                end do
                end if
        end subroutine

end module andersen
