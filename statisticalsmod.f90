module statisticalsmod
use stdmod
contains
        ! Realizado -> Sergio, Revisado -> Alejandro, Victor !
        ! Esta subrutina calcula la desviación estándar de la velocidad !
        ! Para ello lee un archivo con los xyz de la velocidad y trata  !
        ! los datos. Nótese que se necesita saber el promedio.          !
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

        ! Reciclado de otro programa !
        ! Subrutina para hacer un vector de binning !
        subroutine binsvector(nMC,nmeas,timenergyvec,bins,binvec)
                implicit none
                integer, intent(in) :: nmeas
                integer, intent(in) :: nMC, bins 
                real, dimension(nMC/nmeas), intent(in) :: timenergyvec
                real, dimension((nMC/nmeas)/bins), intent(out) :: binvec
                integer, dimension(2) :: binsplit
                integer :: bini
                integer :: newdimension

                newdimension = (nMC/nmeas)/bins

                do bini=1, newdimension
                        binsplit(1) = ((bini-1)*bins)+1
                        binsplit(2) = bini*bins
                        binvec(bini) = sum(timenergyvec(binsplit(1):binsplit(2)))/(newdimension)
                end do
        end subroutine

        ! Reciclado de otro programa !
        ! Subrutina para realizar binning !
        subroutine Binning(NSteps,Mat,Filename) 
                implicit none
                character(6), intent(in) :: Filename
                integer, intent(in) :: NSteps
                real, dimension(NSteps), intent(out) :: Mat
                integer :: BinsN, newdimension
                real, dimension(:), allocatable :: binvec

                BinsN = 1.
                open(unit=30, file=Filename, form='formatted', status='unknown')
                do while (BinsN < NSteps/2)
                        BinsN = BinsN*2
                        newdimension=NSteps/BinsN
                        allocate(binvec(newdimension))
                        call binsvector(NSteps,1,Mat,BinsN,binvec)
                        write(30,*) sdeviation(newdimension,binvec,Nsteps)
                        deallocate(binvec)
                end do
                close(30)
        end subroutine
end module statisticalsmod
