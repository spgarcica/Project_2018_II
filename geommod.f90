module geommod
contains
        !Subrutina que genera la geometria inicial del sistema que corresponde a una estructura cúbica simple!
        subroutine GenGeom(M,Density,CO_R,Position_mat,N_atoms,aa,L_box,Velocity_mat,fuerza,Cutoff)
                implicit none
                integer, intent(in) :: M
	        real, intent(in) :: Density, CO_R
	        real, dimension(:,:), allocatable, intent(out) :: Position_mat, Velocity_mat,fuerza
                integer, intent(out) :: N_atoms
         	real, intent(out) :: aa, L_box, Cutoff
                integer :: ii, jj, kk, cont

                N_atoms = M**3
                L_box = (float(N_atoms)/Density)**(1.0/3.0)
                aa = L_box/M
                Cutoff = CO_R * L_Box
                allocate(Position_mat(N_atoms,3))
                allocate(Velocity_mat(N_atoms,3))
                allocate(fuerza(N_atoms,3))

                cont = 0
                do ii=1, M
                        do jj=1, M
                                do kk=1, M
                                        cont = cont + 1
                                        Position_mat(cont,1) = aa*(ii-1) 
                                        Position_mat(cont,2) = aa*(jj-1)
                                        Position_mat(cont,3) = aa*(kk-1)
                                end do
                        end do
                end do
        end subroutine
end module geommod
