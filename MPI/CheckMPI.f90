module CheckMPI
        contains
                subroutine Check(M,num_proc,myrank)
                        integer, intent(in) :: M, num_proc, myrank
                        if ((mod(M,num_proc) /= 0) .or. (M < num_proc)) then
                                if (myrank == 0) then
                                        print *, "System dimension isn't compatible with process number"
                                end if
                                STOP
                        end if
                end subroutine
end module CheckMPI
