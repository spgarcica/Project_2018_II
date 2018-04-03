module LJ_forcemod

contains
        !------------------------------------------------------------------------------------------------------!
        !  Subroutine that provides a matrix with the forces that operate in each particle and the potential   ! 
        !  energy.                                                                                             !
        !------------------------------------------------------------------------------------------------------!
        subroutine LJ_force(N_atoms,N_Interactions,myrank,num_proc,LJ_IntMat,Position_mat, &
                        Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
               implicit none
               integer, intent(in) :: N_atoms, myrank, num_proc
               integer(8) :: N_Interactions
               integer, dimension(N_Interactions,2), intent(in) :: LJ_IntMat
               real, dimension(N_atoms,3), intent(in) :: Position_mat
               real, intent(out) :: Pot_En
               real, dimension(N_atoms,3), intent(out) :: Force_Mat
               integer(8) :: iii, jjj, kkk, lll, Assigner, A_From, A_Until
               real :: Cutoff2, Distance_2, Distance_4, Distance_6, Distance_8, Distance_12, Distance_14, Dist_Aux
               real, dimension(3) :: Distance_vec, Force_vec
               real :: Cutoff, C_U, Pressure, L_Intend

               Force_Mat = 0.
               Pot_En = 0.
               Cutoff2 = Cutoff**2
               Pressure = 0
               N_Interactions = (N_Atoms*(N_Atoms-1))/2
               Assigner = N_Interactions/num_proc
               A_From = (myrank*Assigner)+1
               A_Until = A_From + Assigner - 1

               do lll=A_From, A_Until
                       ! Determinig particle couple of interaction lll !
                       call Particle_Couple(N_Atoms,lll,iii,jjj)
                       Distance_vec = Position_mat(iii,:)-Position_mat(jjj,:)
                       Distance_2 = Distance_vec(1) ** 2 + Distance_vec(2) ** 2 + Distance_vec(3) ** 2
                       ! If the distance is too high, minumum image conditions are imposed !
                       if (sqrt(Distance_2) > Cutoff) then
                               do kkk=1, 3
                                       Dist_Aux = Position_mat(iii,kkk) - Position_mat(jjj,kkk)
                                       if (abs(Dist_Aux) > cutoff) then
                                               Distance_vec(kkk) = Distance_vec(kkk)-sign(L_Intend,Dist_Aux)
                                       end if
                               end do
                               Distance_2 = Distance_vec(1) ** 2 + Distance_vec(2) ** 2 + Distance_vec(3) ** 2
                       end if
                       Distance_4 = Distance_2**2
                       Distance_6 = Distance_4*Distance_2
                       Distance_8 = Distance_4**2
                       Distance_12 = Distance_6**2
                       Distance_14 = Distance_8*Distance_6
                       if (Distance_2 < Cutoff2) then
                               Force_vec = (48./Distance_14 - 24./Distance_8) * (Distance_vec/Distance_2)
                               Force_Mat(iii,:) = Force_Mat(iii,:) + Force_vec
                               Force_Mat(jjj,:) = Force_Mat(jjj,:) - Force_vec
                               Pot_En = Pot_En + 4.*(1./Distance_12 - 1./Distance_6) - C_U
                               Pressure = Pressure + abs(Dot_Product(Force_vec,Distance_vec))
                       end if
               end do
        end subroutine
        
        subroutine Particle_Couple(N_Atoms,lll,iii,jjj)
                implicit none
                integer, intent(in) :: N_Atoms, lll
                integer, intent(inout) :: iii, jjj
                integer :: incr, counter
                integer(8) :: Accumulator
                
                iii = 0
                jjj = 0
                Accumulator = 0
                counter = 0
                incr = N_Atoms-1
                
                do while (iii == 0 .and. jjj == 0)
                        Accumulator = Accumulator + incr
                        counter = counter + 1
                        if (lll <= Accumulator) then
                                iii = counter
                                jjj = iii + lll - Accumulator + incr
                        end if
                        incr = incr - 1
                end do
        end subroutine
end module LJ_forcemod
