module LJ_forcemod

contains
        ! Realizado -> Sergio, Revisado -> Víctor !
        ! Subrutina que proporciona una matriz con la fuerza que actua en cada partícula y la energia potencial !
        subroutine LJ_force(N_atoms,Position_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
               implicit none
               integer, intent(in) :: N_atoms
               real, dimension(N_atoms,3), intent(in) :: Position_mat
               real, intent(out) :: Pot_En
               real, dimension(N_atoms,3), intent(out) :: Force_Mat
               integer :: iii, jjj, kkk
               real :: Cutoff2, Distance_2, Distance_4, Distance_6, Distance_8, Distance_12, Distance_14, Dist_Aux
               real, dimension(3) :: Distance_vec, Force_vec
               real :: Cutoff, C_U, Pressure, L_Intend

               Force_Mat = 0.
               Pot_En = 0.
               Cutoff2 = Cutoff**2
               Pressure = 0

               do iii=1, N_Atoms
                       do jjj=iii+1, N_Atoms
                               Distance_vec = Position_mat(iii,:)-Position_mat(jjj,:)
                               Distance_2 = Distance_vec(1) ** 2 + Distance_vec(2) ** 2 + Distance_vec(3) ** 2
                               ! Si la distancia es demasiado grande imponemos condiciones de imagen mínima !
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
               end do
        end subroutine
end module LJ_forcemod
