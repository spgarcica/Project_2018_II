module MPI_Regroup
        include 'mpif.h'
        contains
                subroutine Regroup(ierror,myrank,num_proc,N_Slice,N_Dimension,Matrix,R_Matrix)
                        integer, intent(in) :: ierror, myrank, num_proc, N_Slice, N_Dimension
                        real, dimension(N_Slice,N_Dimension), intent(in) :: Matrix
                        real, dimension(N_Slice,N_Dimension), intent(out) :: R_Matrix
                        real, dimension(N_Slice*N_Dimension) :: Aux_Matrix
                        integer :: rank, Counter, ii, jj

                        do rank=0, num_proc-1
                                if (myrank == rank) then
                                        Aux_Matrix = reshape(Matrix,(/N_Dimension*N_Slice/))
                                end if
                                call MPI_BCAST(Aux_Matrix,N_Dimension*N_Slice,MPI_REAL,rank,MPI_COMM_WORLD,ierror)
                                Counter = 1
                                do ii=(rank*N_Slice)+1, (rank+1)*N_Slice
                                        do jj=1, N_Dimension
                                                R_Matrix(ii,jj) = Aux_Matrix(Counter+(N_Slice*(jj-1)))
                                        end do
                                        Counter = Counter + 1
                                end do
                        end do
                end subroutine

                subroutine All_Reduce(ierror,myrank,num_proc,N_Slice,N_Dimension,Matrix,R_Matrix)
                        integer, intent(in) :: ierror, myrank, num_proc, N_Slice, N_Dimension
                        real, dimension(N_Slice*num_proc,N_Dimension), intent(in) :: Matrix
                        real, dimension(N_Slice*num_proc,N_Dimension), intent(out) :: R_Matrix
                        real, dimension(N_Slice*num_proc*N_Dimension) :: Aux_Matrix, Aux_R_Matrix
                        integer :: rank, Counter, ii, jj

                        do rank=0, num_proc-1
                                if (myrank == rank) then
                                        Aux_Matrix = reshape(Matrix,(/N_Dimension*N_Slice*num_proc/))
                                end if
                                call MPI_ALLREDUCE(Aux_Matrix,Aux_R_Matrix,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
                                Counter = 1
                                do ii=1, N_Slice*num_proc
                                        do jj=1, N_Dimension
                                                R_Matrix(ii,jj) = Aux_R_Matrix(Counter+(N_Slice*(jj-1)))
                                        end do
                                        Counter = Counter + 1
                                end do
                        end do
                end subroutine
end module MPI_Regroup
