module repairmod
contains
        subroutine Repair_Cutoff(Cutoff,C_U)
                implicit none
                real, intent(in) :: Cutoff
                real, intent(out) :: C_U
                !He puesto 12 en vez de 2!
                C_U = 4 * (1./Cutoff**12 - 1./Cutoff**6)
        end subroutine
end module repairmod
