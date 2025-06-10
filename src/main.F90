program Isaac

   use iso_c_binding

   !use mpi
   use parameters
   use global
   use parallel_omp
   !use parallel_mpi
   use time
   use output
   use output_data
   use bodies
   use energy
   use quads
   use barnes_hut
   use ini_file_reader
   use input
   use initialize
   use integrators

   implicit none

   interface
        subroutine start_viewer() bind(C)
        end subroutine

        subroutine draw_points(x, y, n) bind(C)
            use iso_c_binding
            real(C_DOUBLE), dimension(*), intent(in) :: x, y
            integer(C_INT), intent(in) :: n
        end subroutine
    end interface

   interface
      subroutine stopFileHandler
      end subroutine stopFileHandler
      
      subroutine getParameterChange
      end subroutine getParameterChange
   end interface
   
   !abstract interface
   !   subroutine evolver ()
   !   end subroutine evolver
   !end interface

   !procedure (evolver), pointer :: evolve => null ()

! ##################################
! VARIABLES
! ##################################

   logical           :: stopFileExists

! ##################################
! START THE PROGRAM
! ##################################

   write (*,*) "############################################################"
   write (*,*) "#   +                   *           *         *  +         #"
   write (*,*) "#   III  *    SSS        AAA     +     AAA           CCC   #"
   write (*,*) "#   III      SSS   *    AA AA         AA AA      +  CCC    #"
   write (*,*) "#   III     SSS        AA   AA   *   AA   AA       CCC     #"
   write (*,*) "#   III      SSS +    AAAAAAAAA     AAAAAAAAA   * CCC   +  #"
   write (*,*) "#   III    +  SSS     AA     AA   + AA     AA      CCC     #"
   write (*,*) "#   III      SSS    * AA   * AA     AA     AA   *   CCC    #"
   write (*,*) "#   III     SSS       AA     AA     AA     AA        CCC   #"
   write (*,*) "#      +          *                   +           +        #"
   write (*,*) "############################################################"
   write (*,*) 
   write (*,*) " ISAAC - barnes-hut gravity simulator, Version 0.4.0 beta"
   write (*,*) " OpenMP parallel version!"
   write (*,*) 
   write (*,*) "------------------------------------------------------------------"

   call start_viewer()

! ##################################
! (A) LOAD THE INPUT DATA
! ##################################

   ! parse the command line for the input file and additional parameters
   call parseCommandLine()
   
   ! check if an input file is given
   if (len_trim(inputFile) .eq. 0) then
      write (*,*) " Missing input file parameter! STOP"
      stop
   end if
   
   ! read the input file
   write (*,'(A,A,A)', advance='no') " Read input file ", trim(inputFile), "... "
   call readInputFile(trim(inputFile)) 
   write (*,*) " done!"
   write (*,*) "------------------------------------------------------------------"

! ##################################
! (B) INITIALIZE THE RUN
! ##################################

   startTimeInit = omp_get_wtime()
      
   write (*,*) " Initialize input parameter variables ..."
   call initializeInputParameters()
   
   write (*,'(A,I2,A)') " Set up ", ompThreads, " OpenMP-Threads ..."
   call omp_set_num_threads(ompThreads)
   
   write (*,*) " Initialize variables ..."
   call initializeVariables()
 
   write (*,*) " Initialize data arrays ..."
   call initializeDataArrays()

   write (*,*) " Initialize bodies ..."
   call initializeBodies()

   write (*,*) " Initialize output ..."
   call initializeOutput()

   ! calculate inital accelerations
   write (*,*) " Calculate initial accelerations ..."
   minimumCollisionTimeSqr_local = calculateAccelerations()
   minimumCollisionTimeSqr       = minimumCollisionTimeSqr_local
   
   ! calculate the initial energies
   write (*,*) " Calculate initial energies ..."
   call calculateEnergies()
   e_in = e_kin + e_pot

   write (*,*) " Finish initialization ..."
   
   if (timestepMode .eq. Timestep_Variable) then
      dt = sqrt(minimumCollisionTimeSqr) * dt_factor
   end if

   ! output for t=0 
   call writeStateFile()
   if (outputEnergyConservation) then
      call writeEnergyFile()
   end if
      
   ! clear up the ini file reader
   call finishIniFileReader()

   ! rename possible stop file
   call rename("STOP", "NOSTOP")

   endTimeInit = omp_get_wtime()

   write (*,*) " Initialization done in ", endTimeInit - startTimeInit," seconds!"
   write (*,*) "------------------------------------------------------------------"

   
! ##################################
! (C) MAIN LOOP
! ##################################

   write (*,*) " Starting the calculation..."
   startTimeCalc = omp_get_wtime()

   ! rev up the multistep integrators
   if ((integrator .eq. Integrator_MultiStep4) .or. (integrator .eq. Integrator_MultiStep6) .or. &
       (integrator .eq. Integrator_MultiStep8)) then
      
      write (*,*) "Rev up multistep integrator ..."
      call rev_up_multistep()
   end if

   ! #####
   ! LOOP
   
   !$omp parallel shared(minimumCollisionTimeSqr) private(minimumCollisionTimeSqr_local)
   do while (.not. reachedEndTime)

      !$omp master
      
      iteration = iteration + 1_16

      dt_old   = dt  
      if (timestepMode .eq. Timestep_Variable) then
         dt = sqrt(minimumCollisionTimeSqr) * dt_factor
      end if
         
      if (outputExactTimes) then
           if (t+dt > t_output) dt = t_output - t
           if (t+dt > t_end)    dt = t_end    - t
      end if
      dt_ratio = dt / dt_old
         
      if (dt .gt. dt_max) dt_max = dt
      if (dt .lt. dt_min) dt_min = dt
         
      if (outputSilentMode) then
         percentage = floor(t/t_end * 100.0d0)
         if (percentage .ne. lastPercentageWritten) then
            write (*,"(I2,A)", advance='no') percentage, "%  "
            lastPercentageWritten = percentage
            if (mod(percentage,10) .eq. 0) write(*,*)
         end if
      else
         write (*, 1000) "iter: ", iteration, "time: ", t+dt, "dt: ", dt, "calc: ", floor(t/t_end * 100.0d0), &
                         "%       dt_ratio: ", dt_ratio
      end if
      
      ! reset the  minimumCollisionTime so we can find the new minimum
      minimumCollisionTimeSqr = huge(0.0d0)
      
      !$omp end master


      !$omp barrier
      
      ! EVOLVE
      select case (integrator)
         case (Integrator_EulerExplicit)
            minimumCollisionTimeSqr_local = evolve_euler_explicit()
         case (Integrator_Leapfrog)
            minimumCollisionTimeSqr_local = evolve_leapfrog()
         case (Integrator_RungeKutta2)
            minimumCollisionTimeSqr_local = evolve_runge_kutta_2()
         case (Integrator_RungeKutta4)
            minimumCollisionTimeSqr_local = evolve_runge_kutta_4()
         case (Integrator_Yoshida4)
            minimumCollisionTimeSqr_local = evolve_yoshida_4()
         case (Integrator_Yoshida6)
            minimumCollisionTimeSqr_local = evolve_yoshida_6()
         case (Integrator_Yoshida8)
            minimumCollisionTimeSqr_local = evolve_yoshida_8()
         case (Integrator_MultiStep4)
            !call evolve_multistep_4()
         case (Integrator_MultiStep6)
            !call evolve_multistep_6()
         case (Integrator_MultiStep8)
            !call evolve_multistep_8()
      end select
      
      ! get the minimum collision time from all threads
      !$omp critical
         minimumCollisionTimeSqr = min(minimumCollisionTimeSqr, minimumCollisionTimeSqr_local)
      !$omp end critical
      
      !$omp master
      call draw_points(body_x, body_y, nrOfBodies)
      !$omp end master

      !$omp master
      t = t + dt
         
      if (t .ge. t_output) then
         
         call calculateEnergies()

         call writeStateFile()

         if (.not. outputSilentMode) then
            write (*,*) "------------------------------------------------------------------"
            write (*,*) " Output for time ", t
            write (*,*) " Relative energy error (e_current - e_in) / e_in = ", (e_kin + e_pot - e_in) / e_in
            write (*,'(A, F6.1)') " Time taken so far: ", omp_get_wtime() - startTimeCalc
            write (*,*) "------------------------------------------------------------------"
         end if

         if (outputEnergyConservation) then
            call writeEnergyFile()
         end if
         
         t_output = t_output + dt_output
      end if
            
      if (outputTimesteps) then
         call writeTimestepFile()
      end if
         
      if (t .ge. t_end) then
         reachedEndTime = .true.
         
         if (outputSilentMode) then
            write (*,'(A)') "100%"
         end if
      end if
         
      ! TODO: rework this to work with MPI
      inquire(file="STOP", exist=stopFileExists)
      
      if (stopFileExists) then
         startTimeUser = omp_get_wtime()
         call stopFileHandler()
         startTimeCalc = startTimeCalc + (omp_get_wtime() - startTimeUser)
      end if   
      !$omp end master
      
      !$omp barrier       
         
   end do
   !$omp end parallel
   
   ! LOOP
   ! #####
   
   endTimeCalc = omp_get_wtime()

! ##################################
! (D) FINAL OUTPUT
! ##################################

   write (*,*) "------------------------------------------------------------------"
   write (*,*) " CALCULATION DONE!"
   write (*,*) "------------------------------------------------------------------"
   write (*,*)

   ! calculate the final energies
   write (*,*) " Calculate final energies ..."
   call calculateEnergies()
   e_out = e_kin + e_pot

   write (*,*) " Write final output ..."
   call writeStateFile()
   call finishOutput()
   
   write (*,*) "------------------------------------------------------------------"
   write (*,*)
   write (*,*) "Energy conservation:"
   write (*,*) "Starting energy:                        e_in = ",          e_in
   write (*,*) "Final energy:                          e_out = ",  e_out
   write (*,*) "Absolute energy error:          e_out - e_in = ",  e_out - e_in
   write (*,*) "Relative energy error: (e_out - e_in) / e_in = ", (e_out - e_in) / e_in
      
   write (*,*)
   write (*,'(A,F16.0)') " Direct force calculations : ", directForceCalculations
   write (*,'(A,F16.0)') " BH Tree force calculations: ", BHTreeForceCalculations
      
   write (*,*)
   write (*,*) "Total # of iterations  = ", iteration
   write (*,*) "Minimum time step used = ", dt_min
   write (*,*) "Average time step used = ", t / dble(iteration)
   write (*,*) "Maximum time step used = ", dt_max
   write (*,*)
      
   write (*,*) "Time for initialization  : ",  endTimeInit - startTimeInit
   write (*,*) "Time for calculation     : ",  endTimeCalc - startTimeCalc
   write (*,*) "Average time per timestep: ", (endTimeCalc - startTimeCalc) / dble(iteration)
   write (*,*) "Total time taken         : ",  endTimeCalc - startTimeInit

   write (*,*) "------------------------------------------------------------------"

! ##################################
! (E) FINALIZE
! ##################################
   
   write (*,*) " Clean up ..."
   
   call resetBHTree()

   deallocate(body_x   ); deallocate(body_y )
   deallocate(body_vx  ); deallocate(body_vy)
   deallocate(body_ax  ); deallocate(body_ay)
   deallocate(body_mass); deallocate(body_mass_inv)
   deallocate(body_mu  )
      
   if (allocated(body_jx )) deallocate(body_jx )
   if (allocated(body_jy )) deallocate(body_jy )
   if (allocated(body_axp)) deallocate(body_axp)
   if (allocated(body_ayp)) deallocate(body_ayp)
   
   write (*,*) " DONE!"  

! ##################################
! (E) FORMATTING STATEMENTS
! ##################################

1000 format(1x, A7, I8, A14, E12.5, A12, F13.3, A11, I2, A18, F6.3)

end program Isaac



! ##################################
! HANDLER IF STOP FILE PRESENT
! ##################################
subroutine stopFileHandler
   use time
   use output_data
   use initialize
   use barnes_hut
   implicit none
   
   integer :: userInput
   
   write (*,*) "STOP file found! Calculation is paused."
   write (*,*) "(1) resume"
   write (*,*) "(2) reload ini file"
   write (*,*) "(3) change parameters"
   write (*,*) "(4) stop calculation"
   write (*,*) "What do you want to do:"
   
   read  (*,'(I1)') userInput
            
   if (userInput .eq. 1) then
      write (*,*) "Renaming the STOP file to NOSTOP and continuing the calculation"
      call rename("STOP", "NOSTOP")
   elseif (userInput .eq. 2) then
      write (*,*) "Ini reload not implemented."
      call rename("STOP", "NOSTOP")
   elseif (userInput .eq. 3) then
      call getParameterChange()
      call rename("STOP", "NOSTOP")
      passParameterChange = .true.
   elseif (userInput .eq. 4) then
      write (*,*) "Renaming the STOP file to NOSTOP and stopping the calculation"
      call rename("STOP", "NOSTOP")
      reachedEndTime = .true.
   else
      write (*,*) "Unkown input. Calculation resumes but 'STOP' file is left present."
   end if
   
end subroutine stopFileHandler


! ##################################
! GET THE PARAMETER CHANGE
! ##################################
subroutine getParameterChange
   use time
   use output_data
   use initialize
   use barnes_hut
   implicit none
   
   integer :: userInput
   
   write (*,*) "Which parameter do you want to change?"
   write (*,*) "(1) t_end"
   write (*,*) "(2) integrator"
   write (*,*) "(3) dt_factor"
   write (*,*) "(4) dt_output"
   write (*,*) "(5) theta"
   write (*,'(A)', advance='no') "Input: "
   read  (*,'(I1)') userInput
   
   select case (userInput)
      case (1)
         write (*,'(A)', advance='no') "t_end: "
         read (*,*) t_end
      case (2)
         write (*,'(A)', advance='no') "integrator: "
         read (*,*) integratorString
         call selectIntegrator()
      case (3)
         write (*,'(A)', advance='no') "dt_factor: "
         read (*,*) dt_factor
      case (4)
         write (*,'(A)', advance='no') "dt_output: "
         read (*,*) dt_output
      case (5)
         write (*,'(A)', advance='no') "theta: "
         read (*,*) theta
         thetaSqr = theta * theta
      case default
         write (*,*) "Unkown input."
   end select
   
end subroutine getParameterChange
