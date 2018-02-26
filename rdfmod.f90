MODULE rdf
   implicit none

   contains
        !-------------------------------------------------------------------------------!
        !          Subroutine that calculates the radial distribution function          !
        !-------------------------------------------------------------------------------!
        subroutine rdf(PosMat,N,dr,NHist,L,Counter1,Nsteps,Density,HistProfile,HistCount)
                implicit none
                integer :: iA, jA, iP
                integer(8), intent(INOUT) :: HistCount
                integer, intent(IN) :: N, NHist, Counter1, Nsteps
                real, intent(IN) :: dR, L, Density
                real :: R, Vol, Current_R, Next_R
                real, parameter :: pi = 3.14159265359
                real, dimension(N,3), intent(IN) :: PosMat
                real, dimension(3) :: dR_vec
                real, dimension(NHist),intent(INOUT) :: HistProfile

                do iA = 1,N-1
                        do jA = iA+1,N
                                dr_vec(1) = PosMat(iA,1) - PosMat(jA,1)
                                dr_vec(2) = PosMat(iA,2) - PosMat(jA,2)
                                dr_vec(3) = PosMat(iA,3) - PosMat(jA,3)
                                call pbc(dR_vec,1,L_intend)

                                R = sqrt(dR_vec(1)*dR_vec(1)+dR_vec(2)*dR_vec(2)+dR_vec(3)*dR_vec(3))
                                if (R .lt. L) then
                                        iP = int(r/dr) + 1
                                        HistProfile(iP) = HistProfile(iP) + 1.0
                                end if
                        end do
                end do

                if(Counter1 .eq. Nsteps) then
                        open(2,file='g_function.txt')
                        do iP = 1,NHist
                                Current_R = (iP-1)*dR
                                Next_R = Current_R + dR
                                Vol = (4.0d0/3.0d0)*pi*((Next_R**3)-(Current_R**3))
                                HistProfile(iP) = (2*HistProfile(iP))/(Vol*Density*N*Nsteps)
                                write(2,*) (iP-0.5)*dR, HistProfile(iP)
                       end do
                       close(2)
                end if
        end subroutine
END MODULE rdf

