module geommod
contains
        !----------------------------------------------------------------------------------------------------!
        !  Subroutine that generates a crystalline structure (simple cubic structure) as an initial geometry ! 
        !----------------------------------------------------------------------------------------------------!
        subroutine GenGeom(M,Density,CO_R,Position_mat,N_atoms,aa,L_box,rank,num_proc,Cutoff)
                implicit none
                integer, intent(in) :: M, rank, num_proc
                integer, intent(out) :: N_atoms
                real, intent(in) :: Density, CO_R
                real, intent(out) :: aa, L_box, Cutoff
                real, dimension(:,:), allocatable, intent(inout) :: Position_mat
                integer :: ii, jj, kk, countador, partition

                N_atoms = M**3
                L_box = (float(N_atoms)/Density)**(1.0/3.0)
                aa = L_box/M
                allocate(Position_mat(N_atoms/num_proc,3))
                Cutoff = CO_R * L_box
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
