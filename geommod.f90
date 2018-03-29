module geommod
contains
        subroutine GenGeom(M,N_atoms,Position_mat,aa,rank,num_proc,aux_vec)
                implicit none
                integer, intent(in) :: M, N_atoms,rank,num_proc
                real, intent(in) :: aa
                real, dimension(:,:), allocatable, intent(inout) :: Position_mat
                real, dimension(:), allocatable, intent(inout) :: aux_vec
                integer :: ii, jj, kk, countador, partition

                allocate(Position_mat(N_atoms/num_proc,3))
                allocate(aux_vec(3*(N_atoms/num_proc)))
                partition=M/num_proc
                countador = 0
                do ii=(rank*partition)+1, (rank+1)*partition
                        do jj=1, M
                                do kk=1, M
                                        countador = countador + 1
                                        Position_mat(countador,1) = aa*(ii-1)
                                        Position_mat(countador,2) = aa*(jj-1)
                                        Position_mat(countador,3) = aa*(kk-1)
                                end do
                        end do
                end do

        end subroutine
end module geommod
