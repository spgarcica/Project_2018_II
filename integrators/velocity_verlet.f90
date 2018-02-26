module velocity_verletmod
   use pbcmod
   use LJ_forcemod
   contains
        ! Realizado -> Cristian !
        ! Algoritmo para Velocity Verlet !
        subroutine velocity_verlet(N,dt,PosMat,VelMat,ForceMat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
                implicit none
                integer :: ii
                integer, intent(in) :: N
                real, intent(in) :: dt
                real, intent(out) :: Pot_En
                real, dimension(N,3) :: NewForceMat
                real, dimension(N,3), intent(inout) :: PosMat, VelMat, ForceMat
                real, intent(in) :: Cutoff, C_U, Pressure, L_Intend

                do ii = 1, N
                        PosMat(ii,:) = PosMat(ii,:) + VelMat(ii,:)*dt + 0.5*ForceMat(ii,:)*dt**2
                end do

                call pbc(PosMat,N,L_Intend)
                call LJ_force(N,PosMat,NewForceMat,Pot_En,Cutoff,C_U,Pressure,L_Intend)

                do ii = 1, N
                        VelMat(ii,:) = VelMat(ii,:) + 0.5*(ForceMat(ii,:) + NewForceMat(ii,:))*dt
                        ForceMat(ii,:) = NewForceMat(ii,:)
                end do
        end subroutine
end module velocity_verletmod
