        ! Realizado -> Cristian !
        ! Algoritmo para Velocity Verlet !
        subroutine velocity_verlet(N,dt,PosMat,VelMat,ForceMat,Pot_En)
                implicit none
                integer :: ii
                integer, intent(in) :: N
                real, intent(in) :: dt
                real, intent(out) :: Pot_En
                real, dimension(N,3) :: NewForceMat
                real, dimension(N,3), intent(inout) :: PosMat, VelMat, ForceMat

                do ii = 1, N
                        PosMat(ii,:) = PosMat(ii,:) + VelMat(ii,:)*dt + 0.5*ForceMat(ii,:)*dt**2
                end do

                call pbc(PosMat,N,L_Intend)
                call LJ_force(N_atoms,PosMat,NewForceMat,Pot_En)

                do ii = 1, N
                        VelMat(ii,:) = VelMat(ii,:) + 0.5*(ForceMat(ii,:) + NewForceMat(ii,:))*dt
                        ForceMat(ii,:) = NewForceMat(ii,:)
                end do
        end subroutine
