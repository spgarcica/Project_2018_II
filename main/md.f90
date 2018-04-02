program MMDyn
        use mtmod
        use normaldist
        use pbcmod
        !use rdfmod
        use LJ_forcemod
        use eulermod
        use velocity_verletmod
        use andersenmod
        use gaussmod
        use geommod
        use initial_velocitymod
        use momentumod
        use repairmod
        use statisticalsmod
        use stdmod 
        use CheckMPI
        include 'mpif.h'
        
        ! MPI values !
        integer :: ierror, num_proc, myrank, status(MPI_STATUS_SIZE), N_Slice
        integer :: M, N_atoms, Seed, NSteps, NSta, ii, jj
        integer(8) :: N_Interactions
        integer, dimension(:,:), allocatable :: LJ_IntMat
        real :: Density, Cutoff, L_Box, L_Intend, dt, A_box, Sigma, Temp, Pot_En, Pressure
        real :: E_Constant, S_Constant, M_Constant, C_U, A_Prob, CO_R, dr
        real, dimension(:,:), allocatable :: Position_mat, Velocity_mat, Force_Mat
        ! MPI values !
        real, dimension(:,:), allocatable :: Position_mat_Total, Velocity_mat_Total, Force_Mat_Total
        real, dimension(:,:), allocatable :: NewForce_Mat_Total
        real, dimension(3) :: partial_sum
        real, dimension(2) :: Time
        character(3) :: Integrator

        call MPI_INIT(ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,num_proc,ierror)
        ! Initialization of the clock !
        call cpu_time(Time(1))
        if ( myrank == 0) then
                print *, 'Reading parameters...'
        end if
        ! Reading the parameters file !
        call Read_file('param.dat')
        if ( myrank == 0) then
                print *, 'Generating geometry...'
        end if
        call Check(M,num_proc,myrank)
        ! Initializing the aleatory numbers generator !
        call sgrnd((myrank+1)*Seed)
        ! Generation of the geometry with the given parameters !
        call GenGeom(M,Density,CO_R,Position_mat,N_atoms,A_box,L_Intend,myrank,num_proc,Cutoff)
        ! Define the size to every processor !
        N_Slice = N_atoms/num_proc
        ! Initialization of MPI !
        allocate(Position_mat_Total(N_atoms,3),Velocity_mat_Total(N_atoms,3),Force_Mat_Total(N_atoms,3))
        allocate(NewForce_mat_total(N_atoms,3))
        ! Generation of the function to repair the Cutoff !
        call Repair_Cutoff(Cutoff,C_U)
        if ( myrank == 0) then
                print *, 'Initializing velocities...'
        end if
        ! Generation of random initial velocities !
        call Initial_Velocity(N_atoms,num_proc,Sigma,Velocity_mat,partial_sum,Force_mat)
        ! Barrier !
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        call Regroup(ierror,myrank,num_proc,N_Slice,3,Position_mat,Position_mat_Total)
        call Regroup(ierror,myrank,num_proc,N_Slice,3,Velocity_mat,Velocity_mat_Total)
        ! Initialization of LJ !
        call Interactions_Init(N_Atoms,LJ_IntMat,N_Interactions,myrank)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        call LJ_force(N_atoms,N_Interactions,myrank,num_proc,LJ_IntMat,Position_mat_Total,&
                Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
        !call LJ_force(N_atoms,Position_mat_Total,Force_Mat_Total,Pot_En,Cutoff,C_U,Pressure,L_Intend)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        call All_Reduce(ierror,myrank,num_proc,N_atoms,3,Force_Mat,Force_Mat_Total)
        ! Calling the function that performs the dynamics !
        call Reaper(N_atoms,NSteps,NSta,Position_mat,Position_mat_Total,Velocity_mat,Velocity_mat_Total, &
                        Force_Mat,Force_Mat_Total,Integrator,myrank,num_proc,N_Slice,ierror)
        !Extracting the time taken to perform the dynamics !
        call cpu_time(Time(2))
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        print *, 'CORE NUMBER', myrank, 'REPORTING  ALL WORK DONE IN:', Time(2)-Time(1), 's'
        call MPI_FINALIZE(ierror)
        contains
        !---------------------------------------------------------------------------------------------------------!
        !                                     Subroutine to read the parameters                                   !
        ! Every parameter is given in reduced units and convert at the end of the program the give the results in !
        ! the wanted units.                                                                                       !
        ! If the cutoff ratio is higher than 0.5, the value won't be accepted due to a problem with high cutoffs  !
        !---------------------------------------------------------------------------------------------------------!
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

        !---------------------------------------------------------------------------------------------!
        !                         Subroutine that performs the molecular dynamic                      !
        !---------------------------------------------------------------------------------------------!
        subroutine Reaper(N_atoms,NSteps,NSta,Position_mat,Position_mat_Total,Velocity_mat,Velocity_mat_Total, &
                        Force_Mat,Force_Mat_Total,Integrator,myrank,num_proc,N_Slice,ierror)
                implicit none
                integer, intent(in) :: N_atoms, NSteps, NSta
                ! MPI Values !
                integer, intent(in) :: myrank, num_proc, N_Slice, ierror
                integer(8) :: HistCount
                character(3), intent(in) :: Integrator
                real, dimension(N_atoms/num_proc,3), intent(inout) :: Position_mat, Velocity_mat
                ! MPI values !
                real, dimension(N_atoms,3), intent(inout) :: Position_mat_Total, Velocity_mat_Total, Force_Mat_Total
                real, dimension(N_atoms,3) :: Force_Mat
                integer :: Counter1, Counter2, Counter3, NHist
                real :: Temperature
                real, dimension(NSteps) :: Pressure_Mat, Temperature_Mat, KE_Mat, PE_Mat
                real, dimension(3) :: Vel_Aver, SD_Value
                real, parameter :: Boltz=1.38064852E-23, Avogadro=6.022140857E23
                real, dimension(:), allocatable :: HistProfile
                procedure (), pointer :: Int_Type => null()

                Temperature = 0
                Vel_Aver = 0.

                ! Pointer to the subroutine with the choosen integrator !
                select case(Integrator)
                case ('Eul')
                        Int_Type => Integrate_Euler
                case ('Ver')
                        Int_Type => Integrate_Verlet
                end select

                if ( myrank == 0) then
                        print *, 'Calculating initialization steps...'
                end if
                do Counter1=1, NSta
                        ! First Andersen with fixed 0.1 probability to stabilize the system !
                        call Int_Type
                end do

                ! Resetting the momentum after the stabilization Andersen !
                do Counter1=1, 3
                        Velocity_mat_Total(:,Counter1) = Velocity_mat_Total(:,Counter1) &
                                - sum(Velocity_mat_Total(:,Counter1))/N_atoms
                end do

                ! Opening of the files where the final data will be outputted !
                if (myrank == 0) then
                        open(unit=25, file='traj.xyz', form='formatted', status='unknown')
                        open(unit=26, file='ener.xyz', form='formatted', status='unknown')
                        open(unit=27, file='testvel', form='formatted', status='unknown')
                        open(unit=28, file='temperature.data', form='formatted', status='unknown')
                        write(26,FMT="(A14,A17,A17,A17,A17,A17)") &
                                'Time', 'Total momentum', 'Pressure', 'Potential Enenergy', 'Kin Energy', 'Total Energy'
                end if

                ! Initializing the g(r) !
                NHist = int(L_intend/dr)
                allocate(HistProfile(NHist+1))
                HistProfile = 0.

                if ( myrank == 0) then
                        print *, 'Calculating the dynamic...'
                end if

                do Counter1=1, NSteps
                        call Int_Type

                        if (myrank == 0) then
                                write(25,*) N_atoms
                                write(25,*) ''
                                do Counter2=1, N_atoms
                                        write(25,*) 'He', Position_mat_Total(Counter2,:)*S_constant
                                end do
                        end if

                        !call rdf(Position_mat,N_atoms,dr,NHist,L_intend,Counter1,Nsteps,Density,HistProfile,HistCount)

                        ! Correcting the units and saving the results !
                        Pressure_Mat(Counter1) = (Pressure/(3*L_Intend**3)+Density*Temp)*(E_Constant/(S_Constant**3))
                        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        call MPI_ALLREDUCE(Pressure_Mat(Counter1),Pressure_Mat(Counter1),1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
                        KE_Mat(Counter1) = Kinetic_E(N_atoms/num_proc,Velocity_mat)*E_Constant
                        call MPI_ALLREDUCE(KE_Mat(Counter1),KE_Mat(Counter1),1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
                        PE_Mat(Counter1) = Pot_En*E_Constant
                        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        call MPI_ALLREDUCE(PE_Mat(Counter1),PE_Mat(Counter1),1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
                        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        Temperature_Mat(Counter1) = (2*KE_Mat(Counter1)/(3.*N_atoms))
                        Temperature = Temperature + Temperature_Mat(Counter1)

                        if (myrank == 0) then
                                write(28,*) Temperature_Mat(Counter1)
                                write(26,*) Counter1*sqrt(E_Constant/(M_Constant*S_Constant**2)),&
                                        p_total(N_atoms,Velocity_mat_Total)*sqrt(M_Constant*E_Constant),&
                                        Pressure_Mat(Counter1), Pot_En,KE_Mat(Counter1), Pot_En+KE_Mat(Counter1)
                        
                                do Counter3=1, N_atoms
                                        write(27,*) Velocity_mat(Counter3,:)*sqrt(E_constant/M_constant)
                                        Vel_Aver = Vel_Aver + Velocity_mat(Counter3,:)*sqrt(E_constant/M_constant)
                                end do
                        end if
                end do
                ! Generation of a Gaussian curve !
                call gauss(Temp*(E_constant/Boltz),N_atoms,NSteps)
                
                ! Printing all the results with their standard deviation !
                if (myrank == 0) then
                        print *, 'Average Temperature (Reduced units)'
                        print *,  Temperature/NSteps, '+-', Sdeviation(NSteps,Temperature_Mat,NSteps)
                        print *, 'Average Pressure'
                        print *,  sum(Pressure_Mat)/NSteps, '+-', Sdeviation(NSteps,Pressure_Mat,Nsteps)
                        print *, 'Average Kinetic Energy'
                        print *,  sum(KE_Mat)/NSteps, '+-', Sdeviation(NSteps,KE_Mat,Nsteps)
                        print *, 'Average Potential Energy'
                        print *,  sum(PE_Mat)/NSteps, '+-', Sdeviation(NSteps,PE_Mat,NSteps)
                        print *, 'Average Total Energy'
                        print *,  sum(PE_Mat+KE_Mat)/NSteps, '+-', Sdeviation(NSteps,PE_Mat+KE_Mat,NSteps)

                        close(25)
                        close(26)
                        close(27)
                end if

                if ( myrank == 0) then
                        Vel_Aver = Vel_Aver/(N_atoms*NSteps)
                        print *, 'Calculating SD of velocities...'
                        call SD_Calculator(N_Atoms*NSteps,'testvel',Vel_Aver,SD_Value)
                        print *, 'Velocity average:', Vel_Aver
                        ! The standard deviation is square root of T !
                        print *, 'Velocity Sdeviation:', SD_Value
                end if
        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!                           MPI Subroutines                           !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine Regroup(ierror,myrank,num_proc,N_Slice,N_Dimension,Matrix,R_Matrix)
                integer, intent(in) :: ierror, myrank, num_proc, N_Slice, N_Dimension
                real, dimension(N_Slice,N_Dimension), intent(in) :: Matrix
                real, dimension(N_Slice*num_proc,N_Dimension), intent(out) :: R_Matrix
                real, dimension(N_Slice*N_Dimension) :: Aux_Matrix
                integer :: rank, Counter, ii, jj

                do rank=0, num_proc-1
                        if (myrank == rank) then
                                Aux_Matrix = reshape(Matrix,(/N_Dimension*N_Slice/))
                        end if
                        call MPI_BCAST(Aux_Matrix,N_Dimension*N_Slice,MPI_REAL,rank,MPI_COMM_WORLD,ierror)
                        Counter = 1
                        do ii=(rank*N_Slice)+1, (rank+1)*N_Slice
                                do jj=1, N_Dimension
                                        R_Matrix(ii,jj) = Aux_Matrix(Counter+(N_Slice*(jj-1)))
                                end do
                                Counter = Counter + 1
                        end do
                end do
        end subroutine

        subroutine All_Reduce(ierror,myrank,num_proc,N_atoms,N_Dimension,Matrix,R_Matrix)
                integer, intent(in) :: ierror, myrank, N_Atoms, N_Dimension, num_proc
                real, dimension(N_atoms,N_Dimension), intent(in) :: Matrix
                real, dimension(N_atoms,N_Dimension), intent(out) :: R_Matrix
                real, dimension(N_atoms*N_Dimension) :: Aux_Matrix, Aux_R_Matrix
                integer :: rank, Counter, ii, jj

                Aux_Matrix = reshape(Matrix,(/N_atoms*N_Dimension/))
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call MPI_ALLREDUCE(Aux_Matrix,Aux_R_Matrix,N_atoms*N_Dimension,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
                R_Matrix = reshape(Aux_R_Matrix,(/N_Atoms,3/))
        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!                       Integration Subroutines                       !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine Integrate_Verlet()
                call Velocity_Verlet(N_atoms,myrank,num_proc,dt,Position_mat_total,& 
                        Position_mat,Velocity_mat_total,Force_mat_total,L_Intend,N_Slice)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call Regroup(ierror,myrank,num_proc,N_Slice,3,Position_mat,Position_mat_Total)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call LJ_force(N_atoms,N_Interactions,myrank,num_proc,LJ_IntMat,Position_mat_Total,&
                        Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call All_Reduce(ierror,myrank,num_proc,N_atoms,3,Force_Mat,NewForce_mat_total)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call Velocity_Verlet2(N_atoms,myrank,num_proc,dt,Velocity_mat_total,&
                        Velocity_mat,Force_mat_total,NewForce_mat_total,N_Slice)
                call Andersen(N_atoms,myrank,num_proc,A_Prob,Velocity_mat,Temp)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call Regroup(ierror,myrank,num_proc,N_Slice,3,Velocity_mat,Velocity_mat_Total)
                Force_Mat_Total = NewForce_Mat_Total
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        end subroutine

        subroutine Integrate_Euler()
                call euler(N_atoms,myrank,num_proc,dt,Position_mat_total,Position_mat,Velocity_mat_total,&
                        Velocity_mat,Force_mat_total,L_Intend,N_Slice)
                call Andersen(N_atoms,myrank,num_proc,A_Prob,Velocity_mat,Temp)
                call Regroup(ierror,myrank,num_proc,N_Slice,3,Position_mat,Position_mat_Total)
                call Regroup(ierror,myrank,num_proc,N_Slice,3,Velocity_mat,Velocity_mat_Total)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call LJ_force(N_atoms,N_Interactions,myrank,num_proc,LJ_IntMat,Position_mat_Total,&
                        Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                call All_Reduce(ierror,myrank,num_proc,N_atoms,3,Force_Mat,Force_mat_total)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        end subroutine
end program
