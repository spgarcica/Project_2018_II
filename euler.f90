        ! Realizado -> Alejandro !
        ! Algoritmo para calcular Euler !
        subroutine Euler(N,dt,PosMat,VelMat,ForceMat,Pot_En)
                implicit none
                integer, intent(in) :: N
	        real, intent(in) :: dt
	        real, dimension(N,3), intent(INOUT) :: PosMat, VelMat, ForceMat
	        real, dimension(N,3) :: MatForce
	        real, intent(out) :: Pot_En
                integer :: ii

                call pbc(PosMat,N,L_Intend)
                call LJ_force(N,PosMat,MatForce,Pot_En)

                do ii = 1, N
                        PosMat(ii,:) = PosMat(ii,:) + VelMat(ii,:)*dt + 0.5*MatForce(ii,:)*dt**2
                end do

                do ii = 1, N
                        VelMat(ii,:) = VelMat(ii,:) + MatForce(ii,:)*dt
                end do
        end subroutine
