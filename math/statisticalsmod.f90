module statisticalsmod
use stdmod
contains
        !----------------------------------------------------------------------------!
        !     Subroutine that calculates the standard deviation of the velocities    !
        ! It reads the velocities from a file and has to know the average of the     !
        ! results listed in the file.                                                !
        !----------------------------------------------------------------------------!
        subroutine SD_Calculator(N_Units,Filename,Average,SD_Value)
                implicit none
                integer, intent(in) :: N_Units
                real, dimension(3), intent(in) :: Average
                real, dimension(3), intent(out) :: SD_Value
                character(7) :: Filename
                integer :: ii
                real, dimension(3) :: Aux_Vec

                SD_Value = 0.

                open(unit=25, file=filename, form='formatted', status='unknown')
                do ii=1, N_Units
                        read(25,*) Aux_Vec
                        SD_Value = SD_Value + (Aux_Vec - Average)**2
                end do
                SD_Value = SD_Value * (1./(N_Units-1.))
                SD_Value = sqrt(SD_Value)
        end subroutine
end module statisticalsmod
