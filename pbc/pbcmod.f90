MODULE pbcmod
   contains
        subroutine pbc(Pos_mat,partition,L)
        implicit none
        integer, intent(inout)  :: partition
        real, dimension(partition,3), intent(inout) :: Pos_mat
        real, intent(in) :: L
        integer :: ii, jj

        do ii=1, partition
                do jj=1,3
                        if (Pos_mat(ii,jj) .lt. 0) then
                                  Pos_mat(ii,jj) = Pos_mat(ii,jj) + PBC_Cor(Pos_mat(ii,jj),L) + L
                        else if (Pos_mat(ii,jj) .gt. L) then
                                  Pos_mat(ii,jj) = Pos_mat(ii,jj) - PBC_Cor(Pos_mat(ii,jj),L)
                        end if
                end do
        end do
        end subroutine

        real function PBC_Cor(Val,Siz)
        implicit none
        real, intent(in) :: Val, Siz
        PBC_Cor = abs(aint(Val/Siz))*Siz
        end function PBC_Cor
END MODULE pbcmod
