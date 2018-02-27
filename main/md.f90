program MMDyn
        use mtmod
        use normaldist
        use pbcmod
        use rdfmod
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
        implicit none
        integer :: M, N_atoms, Seed, NSteps, NSta
        real :: Density, Cutoff, L_Box, L_Intend, dt, A_box, Sigma, Temp, Pot_En, Pressure
        real :: E_Constant, S_Constant, M_Constant, C_U, A_Prob, CO_R, dr
        real, dimension(:,:), allocatable :: Position_mat, Velocity_mat, Force_Mat
        real, dimension(2) :: Time
        character(3) :: Integrator

        ! Initialization of the clock !
        call cpu_time(Time(1))
        print *, 'Reading parameters...'
        ! Reading the parameters file !
        call Read_file('param.dat')
        ! Initializing the aleatory numbers generator !
        call sgrnd(Seed)
        print *, 'Generating geometry...'
        ! Generation of the geometry with the given parameters !
        call GenGeom(M,Density,CO_R,Position_mat,N_atoms,A_box,L_Intend,Velocity_mat,Force_Mat,Cutoff)
        ! Generation of the function to repair the Cutoff !
        call Repair_Cutoff(Cutoff,C_U)
        print *, 'Initializing velocities...'
        ! Generation of random initial velocities !
        call Initial_Velocity(N_atoms,Sigma,Velocity_mat)
        ! Generation of the inital matrix of forces !
        call LJ_force(N_atoms,Position_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
        ! Calling the function that performs the dynamics !
        call Reaper(N_atoms,NSteps,NSta,Position_mat,Velocity_mat,Force_Mat,Integrator)
        ! Extracting the time taken to perform the dynamics !
        call cpu_time(Time(2))
        print *, 'All work done in:', Time(2)-Time(1), 's'
 
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
        subroutine Reaper(N_atoms,NSteps,NSta,Position_mat,Velocity_mat,Force_Mat,Integrator)
                implicit none
                integer, intent(in) :: N_atoms, NSteps, NSta
                character(3), intent(in) :: Integrator
                real, dimension(N_atoms,3), intent(inout) :: Position_mat, Velocity_mat, Force_Mat
                integer :: Counter1, Counter2, Counter3, NHist
                real :: Pot_En, Temperature
                real, dimension(NSteps) :: Pressure_Mat, Temperature_Mat, KE_Mat, PE_Mat
                real, dimension(3) :: Vel_Aver, SD_Value
                real, parameter :: Boltz=1.38064852E-23
                real, dimension(:), allocatable :: HistProfile
                procedure (), pointer :: Int_Type => null()

                Temperature = 0
                Vel_Aver = 0.

                ! Pointer to the subroutine with the choosen integrator !
                select case(Integrator)
                case ('Eul')
                        Int_Type => Euler
                case ('Ver')
                        Int_Type => velocity_verlet
                end select

                print *, 'Calculating initialization steps...'
                do Counter1=1, NSta
                        ! Writes a progression bar. Use -O3 option when compiling !
                        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
                           " Percent Complete: ", (real(Counter1)/real(NSta))*100.0, "%"

                        ! First Andersen with fixed 0.1 probability to stabilize the system !
                        call Andersen(N_atoms,0.1,Velocity_mat,Temp)
                        call Int_Type(N_atoms,dt,Position_mat,Velocity_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
                end do
                print *, ''

                ! Resetting the momentum after the stabilization Andersen !
                do Counter1=1, 3
                        Velocity_mat(:,Counter1) = Velocity_mat(:,Counter1) - sum(Velocity_mat(:,Counter1))/N_atoms
                end do

                ! Opening of the files where the final data will be outputted !
                open(unit=25, file='traj.xyz', form='formatted', status='unknown')
                open(unit=26, file='ener.xyz', form='formatted', status='unknown')
                open(unit=27, file='testvel', form='formatted', status='unknown')
                write(26,FMT="(A14,A17,A17,A17,A17,A17)") &
                        'Time', 'Total momentum', 'Pressure', 'Potential Enenergy', 'Kin Energy', 'Total Energy'

                ! Initializing the g(r) !
                NHist = int(L_intend/dr)
                allocate(HistProfile(NHist+1))
                HistProfile = 0.

                print *, 'Calculating the dynamic...'
                do Counter1=1, NSteps
                        ! Writes a progression bar. Use -O3 option when compiling !
                        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
                           " Percent Complete: ", (real(Counter1)/real(NSteps))*100.0, "%"

                        call Andersen(N_atoms,A_Prob,Velocity_mat,Temp)
                        call Int_Type(N_atoms,dt,Position_mat,Velocity_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)

                        write(25,*) N_atoms
                        write(25,*) ''
                        do Counter2=1, N_atoms
                                write(25,*) 'Ar', Position_mat(Counter2,:)*S_constant
                        end do

                        call rdf(Position_mat,N_atoms,dr,NHist,L_intend,Counter1,Nsteps,Density,HistProfile,HistCount)

                        ! Correcting the units and saving the results !
                        Pressure_Mat(Counter1) = (Pressure/(3*L_Intend**3)+Density*Temp)*(E_Constant/(S_Constant**3))
                        KE_Mat(Counter1) = Kinetic_E(N_atoms,Velocity_mat)*E_Constant
                        PE_Mat(Counter1) = Pot_En*E_Constant
                        Temperature_Mat(Counter1) = (2*KE_Mat(Counter1)/(3.*N_atoms))*(E_constant/Boltz)
                        Temperature = Temperature + Temperature_Mat(Counter1)

                        write(26,*) Counter1*sqrt(E_Constant/(M_Constant*S_Constant**2)),&
                                p_total(N_atoms,Velocity_mat)*sqrt(M_Constant*E_Constant),&
                                Pressure_Mat(Counter1), Pot_En,KE_Mat(Counter1), Pot_En+KE_Mat(Counter1)
                        
                        do Counter3=1, N_atoms
                                write(27,*) Velocity_mat(Counter3,:)*sqrt(E_constant/M_constant)
                                Vel_Aver = Vel_Aver + Velocity_mat(Counter3,:)*sqrt(E_constant/M_constant)
                        end do
                end do
                ! Generation of a Gaussian curve !
                call gauss(Temp*(E_constant/Boltz),N_atoms,NSteps)
                print *, ''
                
                ! Printing all the results with their standard deviation !
                print *, 'Average Temperature (Reduced units)'
                print *,  Temperature/NSteps*Boltz/E_Constant, '+-', Sdeviation(NSteps,Temperature_Mat*Boltz/E_Constant,NSteps)
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

                ! Binning the data to improve the standard deviation calculations !
                call Binning(NSteps,Pressure_Mat,'pressu') 
                call Binning(NSteps,Temperature_Mat,'temper') 
                call Binning(NSteps,KE_Mat,'kineti') 
                call Binning(NSteps,PE_Mat,'potent') 

                Vel_Aver = Vel_Aver/(N_atoms*NSteps)
                print *, 'Calculating SD of velocities...'
                call SD_Calculator(N_Atoms*NSteps,'testvel',Vel_Aver,SD_Value)
                print *, 'Velocity average:', Vel_Aver
                ! The standard deviation is square root of T !
                print *, 'Velocity Sdeviation:', SD_Value
        end subroutine
end program

