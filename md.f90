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

        ! Llamamos al reloj para saber el tiempo que tardaremos !
        call cpu_time(Time(1))
        print *, 'Reading parameters...'
        ! Leemos el archivo de parámetros !
        call Read_file('param.dat')
        ! Inicializamos el generador de número aleatorios !
        call sgrnd(Seed)
        print *, 'Generating geometry...'
        ! Generamos la geometría con los parámetros leidos !
        call GenGeom(M,Density,CO_R,Position_mat,N_atoms,A_box,L_Intend,Velocity_mat,Force_Mat,Cutoff)
        ! Generamos una función para reparar el cutoff !
        call Repair_Cutoff(Cutoff,C_U)
        !Ponemos la matriz en el centro, cambio de sistema de referencia útil al aplicar PBC!
        print *, 'Initializing velocities...'
        ! Generamos las velocidades iniciales al azar !
        call Initial_Velocity(N_atoms,Sigma,Velocity_mat)
        ! Generamos la primera matriz de fuerza !
        call LJ_force(N_atoms,Position_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
        ! Llamamos a la función de la dinámica !
        call Reaper(N_atoms,NSteps,NSta,Position_mat,Velocity_mat,Force_Mat,Integrator)
        ! Volvemos a mirar el reloj !
        call cpu_time(Time(2))
        print *, 'All work done in:', Time(2)-Time(1), 's'
 
        contains 
        ! Realizado -> Sergio !
        !Subrutina donde se especifica al programa lo que se desea hacer y se le aportan los parámetros necesarios!
        subroutine Read_file(filename)
                implicit none
                character(9), intent(in) :: filename

                open(unit=25, file=filename, form='formatted', status='old')
                read(25,*)
                read(25,*) M, Density, CO_R, dt, Seed, Sigma, Temp !Se aporta todo en unidades reducidas!
                ! Error si Cutoff muy grande !
                if (0.5 < CO_R) then
                        print *, 'INPUT ERROR: Invalid cutoff ratio, cutoff ratio must be smaller than 0.5'
                        STOP
                end if
                read(25,*) 
                read(25,*) 
                read(25,*) NSteps, NSta, S_Constant, E_Constant, M_Constant, A_Prob !Se corrigen las unidades después !
                read(25,*)
                read(25,*)
                read(25,*) Integrator, dr
                ! Error si integrador no válido !
                if (Integrator /= 'Ver' .and. Integrator /= 'Eul') then
                        print *, 'INPUT ERROR: Integrator must be Ver o Eul'
                        STOP
                end if
                close(25)
        end subroutine

        ! Realizado -> Sergio, Revisado -> Todos !
        ! Esta es la subrutina principal, se encarga de realizar la dinámica con los datos    !
        ! introducidos y de guardar los archivos, también corrije las unidades y calcula g(r) !
        subroutine Reaper(N_atoms,NSteps,NSta,Position_mat,Velocity_mat,Force_Mat,Integrator)
                implicit none
                integer, intent(in) :: N_atoms, NSteps, NSta
                character(3), intent(in) :: Integrator
                real, dimension(N_atoms,3), intent(inout) :: Position_mat, Velocity_mat, Force_Mat
                integer :: Counter1, Counter2, Counter3, NHist
                integer(8) :: HistCount
                real :: Pot_En, Temperature
                real, dimension(NSteps) :: Pressure_Mat, Temperature_Mat, KE_Mat, PE_Mat
                real, dimension(3) :: Vel_Aver, SD_Value
                real, parameter :: Boltz=1.38064852E-23
                real, dimension(:), allocatable :: HistProfile
                procedure (), pointer :: Int_Type => null()

                Temperature = 0
                Vel_Aver = 0.

                ! Apuntamos a las subrutinas según integrador !
                select case(Integrator)
                case ('Eul')
                        Int_Type => Euler
                case ('Ver')
                        Int_Type => velocity_verlet
                end select

                print *, 'Calculating initialization steps...'
                do Counter1=1, NSta
                        ! Esto escribe una barra de progreso, usar -O3 al compilar !
                        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
                           " Percent Complete: ", (real(Counter1)/real(NSta))*100.0, "%"

                        ! Primer Andersen con 0.1 de probabilidad fija para estabilizar el sistema !
                        call Andersen(N_atoms,0.1,Velocity_mat,Temp)
                        call Int_Type(N_atoms,dt,Position_mat,Velocity_mat,Force_Mat,Pot_En,Cutoff,C_U,Pressure,L_Intend)
                end do
                print *, ''

                ! Corregimos el momento a 0 después del primer Andersen de estabilización !
                do Counter1=1, 3
                        Velocity_mat(:,Counter1) = Velocity_mat(:,Counter1) - sum(Velocity_mat(:,Counter1))/N_atoms
                end do

                ! Abrimos los archivos dónde volvamos los datos !
                open(unit=25, file='traj.xyz', form='formatted', status='unknown')
                open(unit=26, file='ener.xyz', form='formatted', status='unknown')
                open(unit=27, file='testvel', form='formatted', status='unknown')
                write(26,FMT="(A14,A17,A17,A17,A17,A17)") &
                        'Time', 'Total momentum', 'Pressure', 'Potential Enenergy', 'Kin Energy', 'Total Energy'

                ! Inizializamos g(r) !
                NHist = int(L_intend/dr)
                HistCount = 0.
                allocate(HistProfile(NHist+1))
                HistProfile = 0.

                print *, 'Calculating the dynamic...'
                do Counter1=1, NSteps
                        ! Esto escribe una barra de progreso, usar -O3 al compilat !
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

                        ! Corregimos unidades y guardamos !
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
                ! Generamos una gaussiana para después comparar con las velocidades !
                call gauss(Temp*(E_constant/Boltz),N_atoms,NSteps)
                print *, ''
                
                ! Imprimimos los datos con su desvest !
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

                ! Realizamos binning para calcular mejor desvest  !
                ! no funciona demasiado bien con valores pequeños !
                ! de NSteps y NSta                                !
                call Binning(NSteps,Pressure_Mat,'pressu') 
                call Binning(NSteps,Temperature_Mat,'temper') 
                call Binning(NSteps,KE_Mat,'kineti') 
                call Binning(NSteps,PE_Mat,'potent') 

                Vel_Aver = Vel_Aver/(N_atoms*NSteps)
                print *, 'Calculating SD of velocities...'
                call SD_Calculator(N_Atoms*NSteps,'testvel',Vel_Aver,SD_Value)
                print *, 'Velocity average:', Vel_Aver
                ! La desviación estándar es raiz de T !
                print *, 'Velocity Sdeviation:', SD_Value
        end subroutine
end program

