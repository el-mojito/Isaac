program Isaac

   use omp_lib
   use parameters
   use global
   use omp_parallel
   use time
   use output
   use output_data
   use bodies
   use energy
   use quads
   use barnes_hut
   use ini_file_reader

   implicit none

   interface
      subroutine evolve_euler_explicit
      end subroutine evolve_euler_explicit
      
      subroutine evolve_leapfrog
      end subroutine evolve_leapfrog
      
      subroutine evolve_runge_kutta_2
      end subroutine evolve_runge_kutta_2
      
      subroutine evolve_runge_kutta_4
      end subroutine evolve_runge_kutta_4
      
      subroutine evolve_yoshida_4
      end subroutine evolve_yoshida_4
      
      subroutine evolve_yoshida_6
      end subroutine evolve_yoshida_6
      
      subroutine evolve_yoshida_8
      end subroutine evolve_yoshida_8
      
      subroutine rev_up_multistep
      end subroutine rev_up_multistep
      
      subroutine evolve_multistep_4
      end subroutine evolve_multistep_4
      
      subroutine evolve_multistep_6
      end subroutine evolve_multistep_6
      
      subroutine evolve_multistep_8
      end subroutine evolve_multistep_8
   end interface
   
   !abstract interface
   !   subroutine evolver ()
   !   end subroutine evolver
   !end interface

   !procedure (evolver), pointer :: evolve => null ()

! ##################################
! VARIABLES
! ##################################

   integer                  :: i, j
   integer                  :: ios
   character(len=50)        :: commandLineArgument 
   double precision         :: f(2)  
   integer                  :: nrOfInits
   logical                  :: isDefault, fileExists
   
   integer                  :: input
   
   character(len=value_len) :: integratorString, forceString, timestepModeString
   
   integer                  :: nrOfObjects, counterStart
   character(len=value_len) :: initString, initType, rotation
   logical                  :: coIsVirtual, objectIsValid
   double precision         :: centerX, centerY, coMass, velX, velY, radius, radius2, minMass, maxMass, velFactorMin, velFactorMax        
   double precision         :: posX, posY, vel, ang, rad, fac1=1.0d0, fac2=1.0d0

   
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
   write (*,*) " ISAAC - barnes-hut gravity simulator, Version 0.2.1 beta"
   write (*,*) 
   write (*,*) "------------------------------------------------------------------"


   
! ##################################
! (A) LOAD THE INPUT DATA
! ##################################

   call get_command_argument(1, commandLineArgument)
   
   if (len_trim(commandLineArgument) .eq. 0) then
      write (*,*) "Missing input file parameter! STOP"
      stop
   end if
     
   ! check & open the input file
   inquire(file=trim(commandLineArgument), exist=fileExists, iostat=ios)

   if (.not. fileExists) then
      write (*,*) "Input file ", trim(commandLineArgument), " not found! STOP"
      stop
   end if   
   
   open (unit=fid_Input, action="read", file=trim(commandLineArgument), iostat=ios)

   if (ios .ne. 0) then
      write (*,*) "Could not open the input file! STOP"
      stop
   end if
   
   ! analyze the input file
   call initIniFileReader(fid_Input, 300)
   call analyzeIniFile(ios)
   
   if (ios .ne. 0) then
      write (*,*) "There was a problem with the input file! STOP"
      call finishIniFileReader()
      close (fid_Input)
      stop
   end if
  
! ##################################
! (B) INITIALIZE THE RUN
! ##################################

   ! general input data
   projectName           = getValue("project_name"          , "project_name_not_set", isDefault)
   gravitationalConstant = getValue("gravitational_constant", 6.6738480d-11         , isDefault)
   ompThreads            = getValue("omp_threads"           , 1                     , isDefault)
   ompSplit              = getValue("omp_split"             , 1                     , isDefault)
   randomSeed            = getValue("random_seed"           , 1                     , isDefault)
   integratorString      = getValue("integrator"            , "leapfrog"            , isDefault)
   forceString           = getValue("force_calculation"     , "direct"              , isDefault)

   select case (trim(integratorString))
      case ("euler explicit", "euler_explicit", "eulerexplicit")
         integrator = Integrator_EulerExplicit
         !evolve => evolve_euler_explicit
      case ("leapfrog")
         integrator = Integrator_Leapfrog
         !evolve => evolve_leapfrog
      case ("runge kutta 2", "runge_kutta_2", "rungekutta2", "rk2")
         integrator = Integrator_RungeKutta2
         !evolve => evolve_runge_kutta_2
      case ("runge kutta 4", "runge_kutta_4", "rungekutta4", "rk4")
         integrator = Integrator_RungeKutta4
         !evolve => evolve_runge_kutta_4
      case ("yoshida 4", "yoshida_4", "yoshida4", "yo4")
         integrator = Integrator_Yoshida4
         !evolve => evolve_yoshida_4
      case ("yoshida 6", "yoshida_6", "yoshida6", "yo6")
         integrator = Integrator_Yoshida6
         !evolve => evolve_yoshida_6
      case ("yoshida 8", "yoshida_8", "yoshida8", "yo8")
         integrator = Integrator_Yoshida8
         !evolve => evolve_yoshida_8
      case ("multistep 4", "multistep_4", "multistep4", "ms4")
         integrator = Integrator_MultiStep4
         !evolve => evolve_multistep_4
      case ("multistep 6", "multistep_6", "multistep6", "ms6")
         integrator = Integrator_MultiStep6
         !evolve => evolve_multistep_6
      case ("multistep 8", "multistep_8", "multistep8", "ms8")
         integrator = Integrator_MultiStep8
         !evolve => evolve_multistep_8
      case default
         write (*,*) "Unkown integrator selected: ", trim(integratorString), "! STOP"
         call finishIniFileReader()
         close (fid_Input)
         stop
   end select
   
   select case (trim(forceString))
      case ("direct", "d")
         forceCalculation = Force_Direct
      case ("barnes hut", "barnes_hut", "bh")
         forceCalculation = Force_BarnesHut
      case default
         write (*,*) "Unkown force calculation selected: ", trim(forceString), "! STOP"
         call finishIniFileReader()
         close (fid_Input)
         stop
   end select
   
   call random_seed()      !TODO: change random seed

   if ((ompThreads .lt. 1) .or. (ompSplit .lt. 1)) then
      write (*,*) "Number for OpenMP-Threads and workload split have to be > 1. STOP"
      call finishIniFileReader()
      close (fid_Input)
      stop
   end if
   
   if (ompThreads .eq. 1) ompSplit = 1

   call omp_set_num_threads(ompThreads)

   ! time, timestep and output data
   t_end                    = getValue("t_end"                     , 1.0d0     , isDefault)
   dt_user                  = getValue("dt"                        , 0.01d0    , isDefault)
   dt_factor                = getValue("dt_factor"                 , 0.01d0    , isDefault)
   dt_output                = getValue("dt_output"                 , 0.1d0     , isDefault)
   timestepModeString       = getValue("timestep_mode"             , "variable", isDefault)
   outputToSingleFile       = getValue("output_to_single_file"     , .false.   , isDefault)
   outputEnergyConservation = getValue("output_energy_conservation", .true.    , isDefault)
   outputTimesteps          = getValue("output_timesteps"          , .true.    , isDefault)
   outputSilentMode         = getValue("output_silent_mode"        , .false.   , isDefault)

   if (trim(timestepModeString) .eq. "constant") then
      timestepMode          = Timestep_Constant
      dt_CalculationEnable  = .false.
      dt_CalculationDisable = .false.
   elseif (trim(timestepModeString) .eq. "variable") then
      timestepMode          = Timestep_Variable
      dt_CalculationEnable  = .true.
      dt_CalculationDisable = .false.
      if ((integrator .eq. Integrator_MultiStep4) .or. (integrator .eq. Integrator_MultiStep6) .or. &
          (integrator .eq. Integrator_MultiStep8)) then
         write (*,*) "Multistep integrator is not yet compatible with variable timestep. STOP"
         call finishIniFileReader()
         close (fid_Input)
         stop
      end if
   else
      write (*,*) "Unknown paramter for timestep mode. STOP"
      call finishIniFileReader()
      close (fid_Input)
      stop
   end if

   ! Barnes Hut data
   theta                    = getValue("theta"                     , 0.2d0  , isDefault)
   thetaSqr                 = theta * theta
   
   ! how to treat runaway objects
   !doRemoveRunaways                = configReader.getValueOfKey<bool>("remove_runaways", false);
   !removeRunawaysIterations        = 1000;								// supress "could be used un-initialized" warning
   !removeRunawaysDistanceFromCom   = 100;								// supress "could be used un-initialized" warning
   !if (doRemoveRunaways) {
   !  removeRunawaysDistanceFromCom = configReader.getValueOfKey<double>("remove_runaways_distance_from_com", 100.0);
   !  removeRunawaysIterations      = configReader.getValueOfKey<int>   ("remove_runaways", 1000);
   !}  

   maxNrOfBodies = getValue("max_nr_of_objects", 1000, isDefault)
   if (isDefault) then
      write (*,*) "No number for max objects specified! STOP"
      call finishIniFileReader()
      close (fid_Input)
      stop
   end if

   allocate(body_x   (maxNrOfBodies)); allocate(body_y       (maxNrOfBodies))
   allocate(body_vx  (maxNrOfBodies)); allocate(body_vy      (maxNrOfBodies))
   allocate(body_ax  (maxNrOfBodies)); allocate(body_ay      (maxNrOfBodies))
   allocate(body_mass(maxNrOfBodies)); allocate(body_mass_inv(maxNrOfBodies))
   allocate(body_mu  (maxNrOfBodies))
   
   ! for the multistep integrators we need to remember the accelerations on previous timesteps
   if (integrator .eq. Integrator_Multistep4) then
      allocate(body_axp(3,maxNrOfBodies))
      allocate(body_ayp(3,maxNrOfBodies))
   elseif (integrator .eq. Integrator_Multistep6) then
      allocate(body_axp(5,maxNrOfBodies))
      allocate(body_ayp(5,maxNrOfBodies))
   elseif (integrator .eq. Integrator_Multistep8) then
      allocate(body_axp(7,maxNrOfBodies))
      allocate(body_ayp(7,maxNrOfBodies))
   end if
   
   nrOfInits = getvalue("nr_of_inits", 1, isDefault)
   if (isDefault) then
      write (*,*) "Number of inits not specified! STOP"
      call finishIniFileReader()
      close (fid_Input)
      stop
   end if
   
   nrOfBodies               = 0
   
   ! init the bodies
   do i = 1, nrOfInits
      write (initString, "(A,I0,A)") "init_", i, "_"
      
      initType = getValue(trim(initString)//"src", "object", isDefault)
      
      if (trim(initType) .eq. "object") then
         nrOfBodies = nrOfBodies + 1
         body_x       (nrOfBodies) = getValue(trim(initString)//"x"   , 0.0d0, isDefault)
         body_y       (nrOfBodies) = getValue(trim(initString)//"y"   , 0.0d0, isDefault)
         body_vx      (nrOfBodies) = getValue(trim(initString)//"vx"  , 0.0d0, isDefault)
         body_vy      (nrOfBodies) = getValue(trim(initString)//"vy"  , 0.0d0, isDefault)
         body_mass    (nrOfBodies) = getValue(trim(initString)//"mass", 1.0d0, isDefault)
         body_mass_inv(nrOfBodies) = 1.0d0 / body_mass(nrOfBodies)
         body_mu      (nrOfBodies) = gravitationalConstant * body_mass(nrOfBodies)
      end if
      
      if (trim(initType) .eq. "system") then
         
         nrOfObjects  = getValue(trim(initString)//"nr_of_objects"           , 10         , isDefault)
         centerX      = getValue(trim(initString)//"center_x"                , 0.0d0      , isDefault)
         centerY      = getValue(trim(initString)//"center_y"                , 0.0d0      , isDefault)
         coMass       = getValue(trim(initString)//"center_object_mass"      , 1.0d30     , isDefault)
         coIsVirtual  = getValue(trim(initString)//"center_object_is_virtual", .false.    , isDefault)
         velX         = getValue(trim(initString)//"vel_x"                   , 0.0d0      , isDefault)
         velY         = getValue(trim(initString)//"vel_y"                   , 0.0d0      , isDefault)
         rotation     = getValue(trim(initString)//"rotation"                , "clockwise", isDefault)
         radius       = getValue(trim(initString)//"radius"                  , 1.0d10     , isDefault)
         radius2      = getValue(trim(initString)//"radius2"                 , 1.0d12     , isDefault)
         minMass      = getValue(trim(initString)//"min_mass"                , 1.0d20     , isDefault)
         maxMass      = getValue(trim(initString)//"max_mass"                , 1.0d25     , isDefault)
         velFactorMin = getValue(trim(initString)//"vel_factor_min"          , 1.0d0      , isDefault)
         velFactorMax = getValue(trim(initString)//"vel_factor_max"          , 1.0d0      , isDefault)
      
         if ((trim(rotation) .ne. "clockwise") .and. (trim(rotation) .ne. "anti-clockwise") .and. &
             (trim(rotation) .ne. "counter-clockwise") .and. (trim(rotation) .ne. "both")) then
            write (*,*) "Unknown parameter for rotation direction for system ", i, ". Using clockwise."
            rotation = "clockwise"
         end if
         
         if (trim(rotation) .eq. "clockwise") then
            fac1 =  1.0d0; fac2 = -1.0d0
         end if
         if ((trim(rotation) .eq. "anti-clockwise") .or. (trim(rotation) .eq. "counter-clockwise")) then
            fac1 = -1.0d0; fac2 =  1.0d0
         end if
         
         counterStart = 1
         
         if (.not. coIsVirtual) then
            nrOfBodies                = nrOfBodies + 1
            body_x       (nrOfBodies) = centerX
            body_y       (nrOfBodies) = centerY
            body_vx      (nrOfBodies) = velX
            body_vy      (nrOfBodies) = velY
            body_mass    (nrOfBodies) = coMass
            body_mass_inv(nrOfBodies) = 1.0d0 / coMass
            body_mu      (nrOfBodies) = gravitationalConstant * body_mass(nrOfBodies)
            counterStart = 2      
         end if
         
         do j = counterStart, nrOfObjects
            objectIsValid = .false.
            
            do while (.not. objectIsValid)
               posX = getRandomNumber(-radius, radius)
               posY = getRandomNumber(-radius, radius)
               
               if (posX .eq. 0) posX = 1.0d-15
               
               rad = sqrt(posX*posX + posY*posY)
               
               if ((rad .ge. radius2) .and. (rad .le. radius)) objectIsValid = .true.
            end do
            
            vel = sqrt(gravitationalConstant * coMass / rad) * getRandomNumber(velFactorMin, velFactorMax)
         
            ang = atan(abs(posY/posX))

            if ((posX .lt. 0) .and. (posY .gt. 0)) ang = -ang +         PI
            if ((posX .lt. 0) .and. (posY .lt. 0)) ang =  ang +         PI
            if ((posX .gt. 0) .and. (posY .lt. 0)) ang = -ang + 2.0d0 * PI
            
            if (trim(rotation) .eq. "both") then
               if (getRandomNumber(0.0d0, 2.0d0) .gt. 1.0d0) then
                  fac1 =  1.0d0; fac2 = -1.0d0
               else
                  fac1 = -1.0d0; fac2 =  1.0d0
               end if
            end if
            
            nrOfBodies                = nrOfBodies + 1
            body_x       (nrOfBodies) = centerX + posX
            body_y       (nrOfBodies) = centerX + posY
            body_vx      (nrOfBodies) = velX + fac1 * vel * sin(ang)
            body_vy      (nrOfBodies) = velY + fac2 * vel * cos(ang)
            body_mass    (nrOfBodies) = getRandomNumber(minMass, maxMass)
            body_mass_inv(nrOfBodies) = 1.0d0 / body_mass(nrOfBodies)
            body_mu      (nrOfBodies) = gravitationalConstant * body_mass(nrOfBodies)
         
         end do

      end if
   end do
   
   
   ! init the first BH tree
   do_dtCalculation = dt_CalculationEnable
   
   call initBHTree()
   call resetBHTree()
   maxDistance = getMaxXYDistanceFromCenter()
   call initBHTreeQuad(2.01d0 * maxDistance)

   do i = 1, nrOfBodies
      call insertBodyInBHTree(BH_Tree, i)
   end do
      
   call calculateBHTreeMassDistribution()
   
   call calculateEnergies()
   e_in = e_kin + e_pot
   
   ! calculate inital accelerations
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do

   ! init some variables
   reachedEndTime = .false.
   iteration      = 0
   dt             = dt_user
   t              = 0.0d0
   t_output       = dt_output
   dt_min         = 1.0d300
   dt_max         = 0.0d0
   lastPercentageWritten = 0
   
   directForceCalculations = 0.0
   BHTreeForceCalculations = 0.0
   
   ! initialize the output, output t=0
   call initOutput()
   call writeStateFile()
   if (outputEnergyConservation) call writeEnergyFile()
   
   ! clear up the ini file reader
   call finishIniFileReader()
   close (fid_Input)
     
   ! rev up the multistep integrators
   if ((integrator .eq. Integrator_MultiStep4) .or. (integrator .eq. Integrator_MultiStep6) .or. &
       (integrator .eq. Integrator_MultiStep8)) then
      call rev_up_multistep()
   end if

   write (*,*) " Initialization done!"
   write (*,*) "------------------------------------------------------------------"


! ##################################
! (C) MAIN LOOP
! ##################################

   write (*,*) " Starting the calculation..."

   !$omp parallel
   do while (.not. reachedEndTime)
      !$omp master
      iteration = iteration + 1

      dt_old   = dt  
      if (timestepMode .eq. Timestep_Variable) then
         dt = sqrt(minimumCollisionTimeSqr) * dt_factor
         if (t+dt > t_output) dt = t_output - t
         if (t+dt > t_end)    dt = t_end    - t
      end if
      dt_ratio = dt_old / dt
      
      if (dt .gt. dt_max) dt_max = dt
      if (dt .lt. dt_min) dt_min = dt

      if (outputSilentMode) then
         percentage = floor(t/t_end * 100.0)
         if (percentage .ne. lastPercentageWritten) then
            write (*,"(I2,A)", advance='no') percentage, "%  "
            lastPercentageWritten = percentage
            if ((mod(percentage,10) .eq. 0) .or. percentage .eq. 99) write(*,*)
         end if
      else
         write (*, 1000) "iter: ", iteration, "time: ", t+dt, "dt: ", dt, "calc: ", floor(t/t_end * 100.0), "%       dt_ratio: ", &
                         dt_ratio
      end if
      !$omp end master

      ! EVOLVE
      !call evolve()
      select case (integrator)
         case (Integrator_EulerExplicit)
            call evolve_euler_explicit()
         case (Integrator_Leapfrog)
            call evolve_leapfrog()
         case (Integrator_RungeKutta2)
            call evolve_runge_kutta_2()
         case (Integrator_RungeKutta4)
            call evolve_runge_kutta_4()
         case (Integrator_Yoshida4)
            call evolve_yoshida_4()
         case (Integrator_Yoshida6)
            call evolve_yoshida_6()
         case (Integrator_Yoshida8)
            call evolve_yoshida_8()
         case (Integrator_MultiStep4)
            call evolve_multistep_4()
         case (Integrator_MultiStep6)
            call evolve_multistep_6()
         case (Integrator_MultiStep8)
            call evolve_multistep_8()
      end select

      !$omp master
      t = t + dt
      
      if (t .ge. t_output) then
         call writeStateFile()
      
         call calculateEnergies()
         
         if (.not. outputSilentMode) then
            write (*,*) "------------------------------------------------------------------"
            write (*,*) " Output for time ", t
            write (*,*) " Relative energy error (e_current - e_in) / e_in = ", (e_kin + e_pot - e_in) / e_in
            write (*,*) "------------------------------------------------------------------"
         end if
         t_output = t_output + dt_output
         
         if (outputEnergyConservation) call writeEnergyFile()
      end if
      
      if (outputTimesteps) call writeTimestepFile()
      
      ! build the new BH tree
      call resetBHTree()
      maxDistance = getMaxXYDistanceFromCenter()
      call initBHTreeQuad(2.01d0 * maxDistance)
   
      do i = 1, nrOfBodies
         call insertBodyInBHTree(BH_Tree, i)
      end do
      
      call calculateBHTreeMassDistribution()
      
      if (t .ge. t_end) reachedEndTime = .true.
      
      inquire(file="STOP", exist=fileExists)
      
      if (fileExists) then
         write (*,*) "STOP file found! Calculation is paused."
         write (*,*) "(1) resume"
         write (*,*) "(2) reload ini file"
         write (*,*) "(3) stop calculation"
         write (*,*) "What do you want to do:"
         read  (*,'(I1)') input
         
         if (input .eq. 1) then
            write (*,*) "Renaming the STOP file to NOSTOP and continuing the calculation"
            call rename("STOP", "NOSTOP")
            continue
         elseif (input .eq. 2) then
         elseif (input .eq. 3) then
            write (*,*) "Renaming the STOP file to NOSTOP and stopping the calculation"
            call rename("STOP", "NOSTOP")
            reachedEndTime = .true.
         else
            write (*,*) "Unkown input. Calculation resumes but 'STOP' file is left present."
         end if
      end if
      !$omp end master
   end do
   !$omp end parallel

! ##################################
! (D) FINAL OUTPUT
! ##################################

   call writeStateFile()

   call finishOutput()

   call calculateEnergies()
   e_out = e_kin + e_pot
   
   write (*,*) "------------------------------------------------------------------"
   write (*,*) " CALCULATION DONE!"
   write (*,*) "------------------------------------------------------------------"
   
   write (*,*)
   write (*,*) "Energy conservation:"
   write (*,*) "Starting energy:                        e_in = ", e_in
   write (*,*) "Final energy:                          e_out = ", e_out
   write (*,*) "Absolute energy error:          e_out - e_in = ", e_out - e_in
   write (*,*) "Relative energy error: (e_out - e_in) / e_in = ", (e_out - e_in) / e_in
   
   write (*,*)
   write (*,'(A,F16.0)') "Direct force calculations : ", directForceCalculations
   write (*,'(A,F16.0)') "BH Tree force calculations: ", BHTreeForceCalculations
   
   write (*,*)
   write (*,*) "Total # of iterations  = ", iteration
   write (*,*) "Minimum time step used = ", dt_min
   write (*,*) "Average time step used = ", t / iteration
   write (*,*) "Maximum time step used = ", dt_max

   write (*,*) "------------------------------------------------------------------"
   

! ##################################
! (E) FINALIZE
! ##################################

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

1000 format(1x, A7, I8, A15, E12.5, A15, F12.3, A15, I2, A1, F15.4)

end program Isaac



! ##################################
! MAIN LOOP BARNES HUT
! ##################################
subroutine main_loop_Barnes_Hut()

   implicit none
       
end subroutine main_loop_Barnes_Hut



! ##################################
! INTEGRATOR EULER EXPLICIT
! ##################################
subroutine evolve_euler_explicit()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   integer :: i
   double precision :: f(2)
   
   minimumCollisionTimeSqr = 1.0d300
   
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do
   
   do i = 1, nrOfBodies
      body_x (i) = body_x (i) + body_vx(i) * dt
      body_y (i) = body_y (i) + body_vy(i) * dt
      
      body_vx(i) = body_vx(i) + body_ax(i) * dt
      body_vy(i) = body_vy(i) + body_ay(i) * dt
   end do
      
end subroutine evolve_euler_explicit



! ##################################
! INTEGRATOR LEAPFROG
! ##################################
subroutine evolve_leapfrog()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   integer :: i
   double precision :: f(2)
   
   minimumCollisionTimeSqr = 1.0d300

   !$omp do
   do i = 1, nrOfBodies
      body_vx(i) = body_vx(i) + 0.5d0 * body_ax(i) * dt
      body_vy(i) = body_vy(i) + 0.5d0 * body_ay(i) * dt
      body_x (i) = body_x (i) +         body_vx(i) * dt
      body_y (i) = body_y (i) +         body_vy(i) * dt
   end do
   !$omp end do
   
   !$omp do
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do
   !$omp end do

   !$omp do
   do i = 1, nrOfBodies
      body_vx(i) = body_vx(i) + 0.5d0 * body_ax(i) * dt
      body_vy(i) = body_vy(i) + 0.5d0 * body_ay(i) * dt
   end do
   !$omp end do
   
end subroutine evolve_leapfrog



! ##################################
! INTEGRATOR RUNGE KUTTA 2
! ##################################
subroutine evolve_runge_kutta_2()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   !integer :: i
   !double precision :: f(2), x_old, y_old, vx_half, vy_half
   
   minimumCollisionTimeSqr = 1.0d300
   
   !do i = 1, nrOfBodies
   !   x_old = body_x(i)
   !   y_old = body_y(i)
   !   
   !   do_dtCalculation = dt_CalculationDisable
   !   
   !   call calculateForce(BH_Tree, i, f)
   !   
   !   body_ax(i) = f(1) * body_mass_inv(i)
   !   body_ay(i) = f(2) * body_mass_inv(i)
   !
   !   vx_half = body_vx(i) + 0.5d0 * body_ax(i) * dt
   !   vy_half = body_vy(i) + 0.5d0 * body_ay(i) * dt
   !   
   !   body_x(i) = body_x(i) + 0.5d0 * body_vx(i) * dt
   !   body_y(i) = body_y(i) + 0.5d0 * body_vy(i) * dt
   !   
   !   do_dtCalculation = dt_CalculationEnable
   !   
   !   call calculateForce(BH_Tree, i, f)
   !   
   !   body_ax(i) = f(1) * body_mass_inv(i)
   !   body_ay(i) = f(2) * body_mass_inv(i)
   !   
   !   body_vx(i) = body_vx(i) + body_ax(i) * dt
   !   body_vy(i) = body_vy(i) + body_ay(i) * dt
   !   
   !   body_x(i) = x_old + vx_half * dt
   !   body_y(i) = y_old + vy_half * dt
   !   
   !end do

end subroutine evolve_runge_kutta_2



! ##################################
! INTEGRATOR RUNGE KUTTA 4
! ##################################
subroutine evolve_runge_kutta_4()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   !integer :: i
   !double precision :: f(2)
   
   minimumCollisionTimeSqr = 1.0d300

end subroutine evolve_runge_kutta_4



! ##################################
! INTEGRATOR YOSHIDA 4
! ##################################
subroutine evolve_yoshida_4()
   use omp_lib
   use parameters
   use bodies
   use time
   use barnes_hut

   implicit none
   
   interface
      subroutine evolve_leapfrog
      end subroutine evolve_leapfrog
   end interface
   
   double precision :: dt_base
   
   minimumCollisionTimeSqr = 1.0d300
   
   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable

   dt =  1.351207191959657d0 * dt_base
   call evolve_leapfrog()
   dt = -1.702414383919315d0 * dt_base
   call evolve_leapfrog()
  
   do_dtCalculation = dt_CalculationEnable
   
   dt =  1.351207191959657d0 * dt_base
   call evolve_leapfrog()
   
   dt =  dt_base
   
end subroutine evolve_yoshida_4



! ##################################
! INTEGRATOR YOSHIDA 6
! ##################################
subroutine evolve_yoshida_6()
   use omp_lib
   use parameters
   use bodies
   use time
   use barnes_hut

   implicit none
   
   interface
      subroutine evolve_leapfrog
      end subroutine evolve_leapfrog
   end interface
   
   double precision :: dt_base
   
   minimumCollisionTimeSqr = 1.0d300
   
   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable
    
   dt =  0.784513610477560d0 * dt_base
   call evolve_leapfrog()
   dt =  0.235573213359357d0 * dt_base
   call evolve_leapfrog()
   dt = -1.17767998417887d0  * dt_base
   call evolve_leapfrog()
   
   dt =  1.31518632068391d0  * dt_base
   call evolve_leapfrog()
   
   dt = -1.17767998417887d0  * dt_base
   call evolve_leapfrog()
   dt =  0.235573213359357d0 * dt_base
   call evolve_leapfrog()
   
   do_dtCalculation = dt_CalculationEnable
   
   dt =  0.784513610477560d0 * dt_base
   call evolve_leapfrog()

   dt =  dt_base

end subroutine evolve_yoshida_6



! ##################################
! INTEGRATOR YOSHIDA 8
! ##################################
subroutine evolve_yoshida_8()
   use omp_lib
   use parameters
   use bodies
   use time
   use barnes_hut

   implicit none
   
   interface
      subroutine evolve_leapfrog
      end subroutine evolve_leapfrog
   end interface
   
   double precision :: dt_base
   
   minimumCollisionTimeSqr = 1.0d300
   
   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable

   dt =  0.104242620869991d1   * dt_base
   call evolve_leapfrog()
   dt =  0.182020630970714d1   * dt_base
   call evolve_leapfrog()
   dt =  0.157739928123617d0   * dt_base
   call evolve_leapfrog()
   dt =  0.244002732616735d1   * dt_base
   call evolve_leapfrog()
   dt = -0.716989419708120d-2  * dt_base
   call evolve_leapfrog()
   dt = -0.244699182370524d1   * dt_base
   call evolve_leapfrog()
   dt = -0.161582374150097d1   * dt_base
   call evolve_leapfrog()

   dt = -0.17808286265894516e1  * dt_base
   call evolve_leapfrog()

   dt = -0.161582374150097d1   * dt_base
   call evolve_leapfrog()
   dt = -0.244699182370524d1   * dt_base
   call evolve_leapfrog()
   dt = -0.716989419708120d-2  * dt_base
   call evolve_leapfrog()
   dt =  0.244002732616735d1   * dt_base
   call evolve_leapfrog()
   dt =  0.157739928123617d0   * dt_base
   call evolve_leapfrog()
   dt =  0.182020630970714d1   * dt_base
   call evolve_leapfrog()
   
   do_dtCalculation = dt_CalculationEnable
   
   dt =  0.104242620869991d1   * dt_base
   call evolve_leapfrog()

   dt =  dt_base

end subroutine evolve_yoshida_8



! ##################################
! REV UP MULTISTEP INTEGRATORS
! ##################################
subroutine rev_up_multistep()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   interface
      subroutine evolve_yoshida_4
      end subroutine evolve_yoshida_4
      
      subroutine evolve_yoshida_6
      end subroutine evolve_yoshida_6
      
      subroutine evolve_yoshida_8
      end subroutine evolve_yoshida_8
   end interface
   
   double precision, allocatable :: body_start_x (:), body_start_y (:)
   double precision, allocatable :: body_start_vx(:), body_start_vy(:)
   
   ! TODO: variable timesteps?
   
   allocate (body_start_x (nrOfBodies))
   allocate (body_start_y (nrOfBodies))
   allocate (body_start_vx(nrOfBodies))
   allocate (body_start_vy(nrOfBodies))
   
   body_start_x  = body_x
   body_start_y  = body_y
   body_start_vx = body_vx
   body_start_vy = body_vy
   
   dt = -dt
   
   if (integrator .eq. Integrator_MultiStep4) then
      call evolve_yoshida_4()
      body_axp(1,:) = body_ax
      body_ayp(1,:) = body_ay
      
      call evolve_yoshida_4()
      body_axp(2,:) = body_ax
      body_ayp(2,:) = body_ay
      
      call evolve_yoshida_4()
      body_axp(3,:) = body_ax
      body_ayp(3,:) = body_ay
      
   elseif (integrator .eq. Integrator_MultiStep6) then
      call evolve_yoshida_6()
      body_axp(1,:) = body_ax
      body_ayp(1,:) = body_ay
      
      call evolve_yoshida_6()
      body_axp(2,:) = body_ax
      body_ayp(2,:) = body_ay
      
      call evolve_yoshida_6()
      body_axp(3,:) = body_ax
      body_ayp(3,:) = body_ay    
      
      call evolve_yoshida_6()
      body_axp(4,:) = body_ax
      body_ayp(4,:) = body_ay
      
      call evolve_yoshida_6()
      body_axp(5,:) = body_ax
      body_ayp(5,:) = body_ay
      
   elseif (integrator .eq. Integrator_MultiStep8) then
      call evolve_yoshida_8()
      body_axp(1,:) = body_ax
      body_ayp(1,:) = body_ay
      
      call evolve_yoshida_8()
      body_axp(2,:) = body_ax
      body_ayp(2,:) = body_ay
      
      call evolve_yoshida_8()
      body_axp(3,:) = body_ax
      body_ayp(3,:) = body_ay    
      
      call evolve_yoshida_8()
      body_axp(4,:) = body_ax
      body_ayp(4,:) = body_ay
      
      call evolve_yoshida_8()
      body_axp(5,:) = body_ax
      body_ayp(5,:) = body_ay
      
      call evolve_yoshida_8()
      body_axp(6,:) = body_ax
      body_ayp(6,:) = body_ay
      
      call evolve_yoshida_8()
      body_axp(7,:) = body_ax
      body_ayp(7,:) = body_ay
   end if
   
   body_x  = body_start_x
   body_y  = body_start_y
   body_vx = body_start_vx
   body_vy = body_start_vy

   deallocate(body_start_x )
   deallocate(body_start_y )
   deallocate(body_start_vx)
   deallocate(body_start_vy)

   dt = -dt

end subroutine rev_up_multistep



! ##################################
! INTEGRATOR MULTISTEP 4
! ##################################
subroutine evolve_multistep_4()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   integer :: i
   double precision :: f(2), jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr = 1.0d300
   
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do
   
   !!!
   
   do i = 1, nrOfBodies
      jdt (1) = (11.0d0 / 6.0d0) * body_ax(i) - 3.0d0 * body_axp(1,i) + 1.5d0 * body_axp(2,i) - body_axp(3,i) / 3.0d0
      jdt (2) = (11.0d0 / 6.0d0) * body_ay(i) - 3.0d0 * body_ayp(1,i) + 1.5d0 * body_ayp(2,i) - body_ayp(3,i) / 3.0d0
      
      sdt2(1) =           2.0d0  * body_ax(i) - 5.0d0 * body_axp(1,i) + 4.0d0 * body_axp(2,i) - body_axp(3,i)
      sdt2(2) =           2.0d0  * body_ay(i) - 5.0d0 * body_ayp(1,i) + 4.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      cdt3(1) =                    body_ax(i) - 3.0d0 * body_axp(1,i) + 3.0d0 * body_axp(2,i) - body_axp(3,i)
      cdt3(2) =                    body_ay(i) - 3.0d0 * body_ayp(1,i) + 3.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      body_x (i) = body_x (i) + body_vx(i) * dt + (0.5d0 * body_ax(i) + jdt (1) / 6.0d0 + sdt2(1) / 24.0d0) * dt * dt
      body_y (i) = body_y (i) + body_vy(i) * dt + (0.5d0 * body_ay(i) + jdt (2) / 6.0d0 + sdt2(2) / 24.0d0) * dt * dt
      
      body_vx(i) = body_vx(i) + body_ax(i) * dt + (0.5d0 *     jdt(1) + sdt2(1) / 6.0d0 + cdt3(1) / 24.0d0) * dt
      body_vy(i) = body_vy(i) + body_ay(i) * dt + (0.5d0 *     jdt(2) + sdt2(2) / 6.0d0 + cdt3(2) / 24.0d0) * dt
      
      !jdt (1) = 11.0d0 * body_ax(i) - 18.0d0 * body_axp(1,i) + 9.0d0 * body_axp(2,i) - 2.0d0 * body_axp(3,i)
      !jdt (2) = 11.0d0 * body_ay(i) - 18.0d0 * body_ayp(1,i) + 9.0d0 * body_ayp(2,i) - 2.0d0 * body_ayp(3,i)
      !
      !sdt2(1) =  2.0d0 * body_ax(i) -  5.0d0 * body_axp(1,i) + 4.0d0 * body_axp(2,i) -         body_axp(3,i)
      !sdt2(2) =  2.0d0 * body_ay(i) -  5.0d0 * body_ayp(1,i) + 4.0d0 * body_ayp(2,i) -         body_ayp(3,i)
      !
      !cdt3(1) =          body_ax(i) -  3.0d0 * body_axp(1,i) + 3.0d0 * body_axp(2,i) -         body_axp(3,i)
      !cdt3(2) =          body_ay(i) -  3.0d0 * body_ayp(1,i) + 3.0d0 * body_ayp(2,i) -         body_ayp(3,i)
      !
      !body_x (i) = body_x (i) + (body_vx(i) + ((18.0d0 * body_ax(i) + jdt (1) + 1.5d0 * sdt2(1)) / 36.0d0) * dt) * dt
      !body_y (i) = body_y (i) + (body_vy(i) + ((18.0d0 * body_ay(i) + jdt (2) + 1.5d0 * sdt2(2)) / 36.0d0) * dt) * dt
      !
      !body_vx(i) = body_vx(i) + ((24.0d0 * body_ax(i) + 2.0d0 * jdt(1) + 4.0d0 * sdt2(1) + cdt3(1)) / 24.0d0) * dt
      !body_vy(i) = body_vy(i) + ((24.0d0 * body_ay(i) + 2.0d0 * jdt(2) + 4.0d0 * sdt2(2) + cdt3(2)) / 24.0d0) * dt
      
      !body_axp(3,i) = body_axp(2,i)
      !body_ayp(3,i) = body_ayp(2,i)
      !
      !body_axp(2,i) = body_axp(1,i)
      !body_ayp(2,i) = body_ayp(1,i)
      !
      !body_axp(1,i) = body_ax (  i)
      !body_ayp(1,i) = body_ay (  i)
   end do
   
   body_axp(3,:) = body_axp(2,:)
   body_ayp(3,:) = body_ayp(2,:)
   
   body_axp(2,:) = body_axp(1,:)
   body_ayp(2,:) = body_ayp(1,:)
   
   body_axp(1,:) = body_ax
   body_ayp(1,:) = body_ay

end subroutine evolve_multistep_4



! ##################################
! INTEGRATOR MULTISTEP 4
! ##################################
subroutine evolve_multistep_6()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   integer :: i
   double precision :: f(2), jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr = 1.0d300
   
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do
   
   do i = 1, nrOfBodies
      jdt (1) = (11.0d0 / 6.0d0) * body_ax(i) - 3.0d0 * body_axp(1,i) + 1.5d0 * body_axp(2,i) - body_axp(3,i) / 3.0d0
      jdt (2) = (11.0d0 / 6.0d0) * body_ay(i) - 3.0d0 * body_ayp(1,i) + 1.5d0 * body_ayp(2,i) - body_ayp(3,i) / 3.0d0
      
      sdt2(1) =           2.0d0  * body_ax(i) - 5.0d0 * body_axp(1,i) + 4.0d0 * body_axp(2,i) - body_axp(3,i)
      sdt2(2) =           2.0d0  * body_ay(i) - 5.0d0 * body_ayp(1,i) + 4.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      cdt3(1) =                    body_ax(i) - 3.0d0 * body_axp(1,i) + 3.0d0 * body_axp(2,i) - body_axp(3,i)
      cdt3(2) =                    body_ay(i) - 3.0d0 * body_ayp(1,i) + 3.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      body_x (i) = body_x (i) + body_vx(i) * dt + (0.5d0 * body_ax(i) + jdt (1) / 6.0d0 + sdt2(1) / 24.0d0) * dt * dt
      body_y (i) = body_y (i) + body_vy(i) * dt + (0.5d0 * body_ay(i) + jdt (2) / 6.0d0 + sdt2(2) / 24.0d0) * dt * dt
      
      body_vx(i) = body_vx(i) + body_ax(i) * dt + (0.5d0 *     jdt(1) + sdt2(1) / 6.0d0 + cdt3(1) / 24.0d0) * dt
      body_vy(i) = body_vy(i) + body_ay(i) * dt + (0.5d0 *     jdt(2) + sdt2(2) / 6.0d0 + cdt3(2) / 24.0d0) * dt
      
      body_axp(3,i) = body_axp(2,i)
      body_ayp(3,i) = body_ayp(2,i)
      
      body_axp(2,i) = body_axp(1,i)
      body_ayp(2,i) = body_ayp(1,i)
      
      body_axp(1,i) = body_ax (  i)
      body_ayp(1,i) = body_ay (  i)
   end do

end subroutine evolve_multistep_6



! ##################################
! INTEGRATOR MULTISTEP 4
! ##################################
subroutine evolve_multistep_8()
   use omp_lib
   use parameters
   use bodies
   use barnes_hut

   implicit none
   
   integer :: i
   double precision :: f(2), jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr = 1.0d300
   
   do i = 1, nrOfBodies
      call calculateForce(BH_Tree, i, f)
      
      body_ax(i) = f(1) * body_mass_inv(i)
      body_ay(i) = f(2) * body_mass_inv(i)
   end do
   
   do i = 1, nrOfBodies
      jdt (1) = (11.0d0 / 6.0d0) * body_ax(i) - 3.0d0 * body_axp(1,i) + 1.5d0 * body_axp(2,i) - body_axp(3,i) / 3.0d0
      jdt (2) = (11.0d0 / 6.0d0) * body_ay(i) - 3.0d0 * body_ayp(1,i) + 1.5d0 * body_ayp(2,i) - body_ayp(3,i) / 3.0d0
      
      sdt2(1) =           2.0d0  * body_ax(i) - 5.0d0 * body_axp(1,i) + 4.0d0 * body_axp(2,i) - body_axp(3,i)
      sdt2(2) =           2.0d0  * body_ay(i) - 5.0d0 * body_ayp(1,i) + 4.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      cdt3(1) =                    body_ax(i) - 3.0d0 * body_axp(1,i) + 3.0d0 * body_axp(2,i) - body_axp(3,i)
      cdt3(2) =                    body_ay(i) - 3.0d0 * body_ayp(1,i) + 3.0d0 * body_ayp(2,i) - body_ayp(3,i)
      
      body_x (i) = body_x (i) + body_vx(i) * dt + (0.5d0 * body_ax(i) + jdt (1) / 6.0d0 + sdt2(1) / 24.0d0) * dt * dt
      body_y (i) = body_y (i) + body_vy(i) * dt + (0.5d0 * body_ay(i) + jdt (2) / 6.0d0 + sdt2(2) / 24.0d0) * dt * dt
      
      body_vx(i) = body_vx(i) + body_ax(i) * dt + (0.5d0 *     jdt(1) + sdt2(1) / 6.0d0 + cdt3(1) / 24.0d0) * dt
      body_vy(i) = body_vy(i) + body_ay(i) * dt + (0.5d0 *     jdt(2) + sdt2(2) / 6.0d0 + cdt3(2) / 24.0d0) * dt
      
      body_axp(3,i) = body_axp(2,i)
      body_ayp(3,i) = body_ayp(2,i)
      
      body_axp(2,i) = body_axp(1,i)
      body_ayp(2,i) = body_ayp(1,i)
      
      body_axp(1,i) = body_ax (  i)
      body_ayp(1,i) = body_ay (  i)
   end do

end subroutine evolve_multistep_8
