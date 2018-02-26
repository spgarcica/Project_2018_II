module stdmod
contains
        ! Reutilizado de otro programa: Calcula la Desviación Estándar !
        real function Sdeviation(Nsize,vec,NSteps)
                 implicit none
                 integer, intent(in) :: Nsize, NSteps
                 real, dimension(Nsize), intent(in) :: vec
                 real :: sdaver
                 integer :: sdcounter

                 sdaver = (sum(vec)/Nsize)
                 sdeviation = 0.

                 do sdcounter=1, Nsize
                         sdeviation = sdeviation + (vec(sdcounter)-sdaver)**2
                 end do
                 sdeviation = sqrt((sdeviation)/(float(NSteps-1)))
        end function sdeviation
end module stdmod
