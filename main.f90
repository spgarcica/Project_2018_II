program inicializacion
	use mtmod
	include 'mpif.h'
    	integer :: M, N_atoms, Seed, NSteps, NSta
    	real :: Density, Cutoff, L_Box, L_Intend, dt, A_box, Sigma, Temp, Pot_En, Pressure, partial_EK, total_EK
    	real :: E_Constant, S_Constant, M_Constant, C_U, A_Prob, CO_R, dr
	real, dimension(2) :: Time
	character(3) :: Integrator
	real, dimension(:,:), allocatable :: Pos_mat, Vel_mat, Pos_mat_total, Vel_mat_total, Force_mat, Force_mat_total
	real, dimension(:), allocatable :: aux_pos, aux_vel, aux_pos_total, aux_vel_total
	real, dimension(3) :: partial_sum, total_sum, partial_P, total_P
	integer :: ierror, rank, num_proc, myrank,contador, i, status(MPI_STATUS_SIZE), param
	
    	! Reading the parameters file !
    	call Read_file('param.dat')
    	! Initializing the aleatory numbers generator !
    	call sgrnd(Seed)
    	N_atoms = M**3
    	L_box = (float(N_atoms)/Density)**(1.0/3.0)
    	A_box = L_box/M
   	Cutoff = CO_R * L_Box

	!A la hora de hacer los cálculos todos los procesos sacaran la información de aquí!
	allocate(Pos_mat_total(N_atoms,3),Vel_mat_total(N_atoms,3),Force_mat_total(N_atoms,3))
	!Matrices auxiliares para pasar información entre procesos!
	allocate(aux_pos_total(3*N_atoms),aux_vel_total(3*N_atoms))	

	call MPI_INIT(ierror)
	call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,num_proc,ierror)


	!Comprobación para que no haya problemas entre las dimensiones del sistema y el número de procesadores!
	if ((mod(M,num_proc) /= 0) .or. (M < num_proc)) then
		if (myrank==0) then
			print *, "System dimension isn't compatible with process number"
		end if
		STOP
	end if
	if (myrank==0) then
		print *, 'L:', L_box
		print *, 'Reading parameters...'
		print *, 'Generating geometry...'
	end if
	!Calculo de un parámetro auxiliar!
	param = N_atoms/num_proc
	! Generation of the geometry with the given parameters !
	call GenGeom(M,N_atoms,Pos_mat,A_box,myrank,num_proc,aux_pos)
	!Calculo de velocidades!
	call Initial_Velocity(N_atoms,num_proc,Sigma,Vel_mat,partial_sum,aux_vel,Force_mat)
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)
	
	!Montamos la matriz en cada proceso!
	do rank=0, num_proc-1
		if (myrank==rank) then
			aux_pos=reshape(Pos_mat,(/3*param/))
			aux_vel=reshape(Vel_mat,(/3*param/))
		end if
		call MPI_BCAST(aux_pos,3*param,MPI_REAL,rank,MPI_COMM_WORLD,ierror)
		call MPI_BCAST(aux_vel,3*param,MPI_REAL,rank,MPI_COMM_WORLD,ierror)
		contador = 1
		do i=(rank*param)+1, (rank+1)*param
			Pos_mat_total(i,1)=aux_pos(contador)
			Pos_mat_total(i,2)=aux_pos(contador+param)
			Pos_mat_total(i,3)=aux_pos(contador+2*param)
			Vel_mat_total(i,1)=aux_vel(contador)
			Vel_mat_total(i,2)=aux_vel(contador+param)
			Vel_mat_total(i,3)=aux_vel(contador+2*param)
			contador=contador+1
		end do
	end do
	call MPI_BARRIER(MPI_COMM_WORLD,ierror)
	!Le quitamos la velocidad al centro de masas!
	do i=1, 3
		call MPI_Reduce(partial_sum(i), total_sum(i), 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
	end do
	call MPI_BCAST(total_sum,3,MPI_REAL,0,MPI_COMM_WORLD,ierror)
	total_sum=total_sum/N_atoms
	do i=1, N_atoms
		Vel_mat_total(i,:) = Vel_mat_total(i,:) - total_sum(:)
	end do

	!Comprovacion de que todo ha ido bien!
	if (myrank==1) then
		print *, 'Pos:'
		do i=1, N_atoms
			print *, Pos_mat_total(i,:)
		end do
		print *, 'Vel:'
		do i=1, N_atoms
			print *, Vel_mat_total(i,:)
		end do
	end if
	
	!Esto es para comprobar que todo funciona correctamente!
	Vel_mat_total=1.0
	Pos_mat_total=1.0
	Force_mat_total=1.0
	if (myrank==1) then
		print *, 'Vel:'
		do i=1, N_atoms
			print *, Vel_mat_total(i,:)
		end do
	end if

	!Como calcular la energia cinética en paralelo!
	partial_EK = Kinetic_E(N_atoms,Vel_mat_total,myrank,param)
	call MPI_BARRIER(MPI_COMM_WORLD,ierror)
	call MPI_Reduce(partial_EK, total_EK, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
	call MPI_BCAST(total_EK,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
	if (myrank==1) then
		print *, 'EK_total', total_EK
	end if
	!Como calcular el momento total del sistema en paralelo!
	call P_sum(N_atoms,Vel_mat_total,myrank,param,partial_P)
	call MPI_BARRIER(MPI_COMM_WORLD,ierror)
	do i=1, 3
		call MPI_Reduce(partial_P(i), total_P(i), 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
	end do
	call MPI_BCAST(total_P,3,MPI_REAL,0,MPI_COMM_WORLD,ierror)	
	if (myrank==1) then
		print *, 'P_total:'
		do i=1, 3
			print *, total_P(i)
		end do
	end if
	
	call Velocity_Verlet(N_atoms,myrank,num_proc,dt,Pos_mat_total,Pos_mat,Vel_mat_total,Vel_mat,Force_mat_total,Force_mat,L_box)
	call MPI_BARRIER(MPI_COMM_WORLD,ierror)
	do rank=0, num_proc-1
		if (myrank==rank) then
			aux_pos=reshape(Pos_mat,(/3*param/))
		end if
		call MPI_BCAST(aux_pos,3*param,MPI_REAL,rank,MPI_COMM_WORLD,ierror)
		contador = 1
		do i=(rank*param)+1, (rank+1)*param
			Pos_mat_total(i,1)=aux_pos(contador)
			Pos_mat_total(i,2)=aux_pos(contador+param)
			Pos_mat_total(i,3)=aux_pos(contador+2*param)
			contador=contador+1
		end do
	end do

	Pos_mat = 11.5
	call pbc(Pos_mat,param,L_box)

        if (myrank==1) then
                print *, 'Pos:'
                do i=1, param
                        print *, Pos_mat(i,:)
                end do
        end if
	
	call MPI_FINALIZE(ierror)
	contains
	
	subroutine Read_file(filename)
    	implicit none
    	character(9), intent(in) :: filename

    	open(unit=25, file=filename, form='formatted', status='old')
        read(25,*)
        read(25,*) M, Density, CO_R, dt, Seed, Sigma, Temp       
		! Warning for cutoffs higher than 0.5 !
        if (0.5 < CO_R) then
        	print *, 'INPUT ERROR: Invalid cutoff ratio, cutoff ratio must be smaller than 0.5'
            STOP
        end if
        read(25,*) 
        read(25,*) 
        read(25,*) NSteps, NSta, S_Constant, E_Constant, M_Constant, A_Prob
        read(25,*)
        read(25,*)
        read(25,*) Integrator, dr
        ! Warning if the integrator given is not valid !
        if (Integrator /= 'Ver' .and. Integrator /= 'Eul') then
        	print *, 'INPUT ERROR: Integrator must be Ver o Eul'
            STOP
        end if
        close(25)
        end subroutine

	subroutine GenGeom(M,N_atoms,Position_mat,aa,rank,num_proc,aux_vec)
		implicit none
		integer, intent(in) :: M, N_atoms,rank,num_proc
		real, intent(in) :: aa
		real, dimension(:,:), allocatable, intent(inout) :: Position_mat
		real, dimension(:), allocatable, intent(inout) :: aux_vec
		integer :: ii, jj, kk, countador, partition

		allocate(Position_mat(N_atoms/num_proc,3))
		allocate(aux_vec(3*(N_atoms/num_proc)))
		partition=M/num_proc
		countador = 0
		do ii=(rank*partition)+1, (rank+1)*partition
			do jj=1, M
				do kk=1, M
					countador = countador + 1
					Position_mat(countador,1) = aa*(ii-1)
					Position_mat(countador,2) = aa*(jj-1)
					Position_mat(countador,3) = aa*(kk-1)
				end do
			end do
		end do

	end subroutine

	subroutine Initial_Velocity(N_atoms,num_proc,Sigma,Velocity_mat,partial_sum,aux_vec,Force_mat)
		implicit none
		integer, intent(in) :: N_atoms, num_proc
		real, intent(in) :: Sigma
		real, dimension(:,:), allocatable, intent(inout) :: Velocity_mat, Force_mat
		real, dimension(:), allocatable, intent(inout) :: aux_vec
		real, dimension(3), intent(out) :: partial_sum
		real, dimension(3) :: random
		integer :: ii

		allocate(Velocity_mat(N_atoms/num_proc,3))
		allocate(Force_mat(N_atoms/num_proc,3))
		allocate(aux_vec(3*(N_atoms/num_proc)))
		Velocity_mat = 1
		do ii=1, N_atoms/num_proc
			random(1) = real(grnd())
			random(2) = real(grnd())
			random(3) = real(grnd())
			Velocity_mat(ii,:) = Velocity_mat(ii,:) + random*Sigma
		end do
		do ii=1, 3
			partial_sum(ii) = sum(Velocity_mat(:,ii))
		end do
	end subroutine	

	real function Kinetic_E(N_atoms,Vel_mat_total,rank,param)
	implicit none
	integer, intent(in) :: N_atoms, rank, param
	real :: Norma
	integer :: ii
	real, dimension(N_atoms,3), intent(in) :: Vel_mat_total	

	Kinetic_E = 0.0
	do ii=(rank*param)+1, (rank+1)*param
		Norma = Vel_mat_total(ii,1)**2 + Vel_mat_total(ii,2)**2 + Vel_mat_total(ii,3)**2
		Kinetic_E = Kinetic_E + Norma
	end do
	Kinetic_E = 0.5*Kinetic_E

	end function

	subroutine P_sum(N_atoms,Vel_mat_total,rank,param,partial_P)
	implicit none
	real, dimension(N_atoms,3), intent(in) :: Vel_mat_total
	real, dimension(3), intent(inout) :: partial_P
	integer, intent(in) :: N_atoms, rank, param
	integer :: ii
	
	partial_P = 0.0
	do ii=(rank*param)+1, (rank+1)*param
		partial_P(1) = partial_P(1) + Vel_mat_total(ii,1)
		partial_P(2) = partial_P(2) + Vel_mat_total(ii,2)
		partial_P(3) = partial_P(3) + Vel_mat_total(ii,3)
	end do
	end subroutine

	subroutine Velocity_Verlet(N_atoms,myrank,num_proc,dt,Pos_mat_total,Pos_mat,Vel_mat_total,Vel_mat,Force_mat_total,Force_mat,L)
	implicit none
                integer :: ii, jj
                integer :: partition
                integer, intent(in) :: N_atoms, myrank, num_proc
                real, intent(in) :: dt, L
                real, dimension(N_atoms,3), intent(in) :: Pos_mat_total, Vel_mat_total, Force_mat_total
                real, dimension(N_atoms/num_proc,3), intent(inout) :: Pos_mat, Vel_mat, Force_mat
		partition = N_atoms/num_proc
		jj=0
		do ii=(partition*myrank)+1, (myrank+1)*partition
			jj = jj + 1
			Pos_mat(jj,:) = Pos_mat_total(ii,:) + Vel_mat_total(ii,:)*dt + 0.5*Force_mat_total(ii,:)*dt**2
		end do
		call pbc(Pos_mat,partition,L)
	end subroutine

	subroutine pbc(Pos_mat,partition,L)
	implicit none
        integer, intent(inout)  :: partition
	real, dimension(partition,3), intent(inout) :: Pos_mat
	real, intent(in) :: L
	integer :: ii, jj
	
	do ii=1, partition
        	do jj=1,3
                	if (Pos_mat(ii,jj) .lt. 0) then
                                  Pos_mat(ii,jj) = Pos_mat(ii,jj) + PBC_Cor(Pos_mat(ii,jj),L) + L
                        else if (Pos_mat(ii,jj) .gt. L) then
                                  Pos_mat(ii,jj) = Pos_mat(ii,jj) - PBC_Cor(Pos_mat(ii,jj),L)
                        end if
                end do
        end do
	end subroutine

	real function PBC_Cor(Val,Siz)
	implicit none
	real, intent(in) :: Val, Siz
	PBC_Cor = abs(aint(Val/Siz))*Siz
	end function PBC_Cor

	
end program
