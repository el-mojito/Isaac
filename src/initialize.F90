! ##################################
! INITIALIZATION ROUTINES
! ##################################
module initialize
   !use mpi
   use parameters
   use global
   !use parallel_mpi
   use parallel_omp
   use time
   use output_data
   use bodies
   use energy
   use barnes_hut
   use output
   use input
   use ini_file_reader

   implicit none

   contains
   
   !------------------------------
   ! Initialize the variables from the input file
   !------------------------------   
   subroutine initializeInputParameters
      implicit none
      
      projectName              = getInputValue("project_name"              , "project_name_not_set")
      gravitationalConstant    = getInputValue("gravitational_constant"    , 6.6738480d-11         )
      
      randomSeed               = getInputValue("random_seed"               , 1                     )
      integratorString         = getInputValue("integrator"                , "leapfrog"            )
      forceString              = getInputValue("force_calculation"         , "direct"              )
      
      if (cl_ompThreads .eq. 0) then
         ompThreads            = getInputValue("omp_threads"               , 1                     )
      else
         ompThreads            = cl_ompThreads
      end if
      
      t_end                    = getInputValue("t_end"                     , 1.0d0     )
      dt_user                  = getInputValue("dt"                        , 0.01d0    )
      dt_factor                = getInputValue("dt_factor"                 , 0.01d0    )
      dt_output                = getInputValue("dt_output"                 , 0.1d0     )
      timestepModeString       = getInputValue("timestep_mode"             , "variable")
      
      outputExactTimes         = getInputValue("output_exact_times"        , .true.    )
      outputToSingleFile       = getInputValue("output_to_single_file"     , .false.   )
      outputEnergyConservation = getInputValue("output_energy_conservation", .true.    )
      outputTimesteps          = getInputValue("output_timesteps"          , .true.    )
      outputSilentMode         = getInputValue("output_silent_mode"        , .false.   )
      
      theta                    = getInputValue("theta"                     , 0.2d0)

      call getNumberOfBodies()
      
   end subroutine initializeInputParameters
      
   !------------------------------
   ! Initialize
   !------------------------------   
   subroutine initializeVariables
      
      implicit none
      
      call random_seed()               !TODO: change random seed
      
      call selectIntegrator()
      
      select case (trim(forceString))
         case ("direct", "d")
            forceCalculation = Force_Direct
         case ("barnes hut", "barnes_hut", "bh")
            forceCalculation = Force_BarnesHut
         case default
            write (*,*) "Unkown force calculation selected: ", trim(forceString), "! STOP"
            call finishIniFileReader()
            stop
      end select
           

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
            !call MPI_FINALIZE (mpiErr) 
            stop
         end if
      else
         write (*,*) "Unknown paramter for timestep mode. STOP"
         call finishIniFileReader()
         !call MPI_FINALIZE (mpiErr) 
         stop
      end if
      
      passParameterChange           = .false.
      thetaSqr                      = theta * theta
      do_dtCalculation              = dt_CalculationEnable
      reachedEndTime                = .false.
      iteration                     = 0_16
      dt                            = dt_user
      t                             = 0.0d0
      t_output                      = dt_output
      dt_min                        = huge(0.0d0)
      dt_max                        = 0.0d0
      lastPercentageWritten         = 0
      minimumCollisionTimeSqr       = huge(0.0d0)
      minimumCollisionTimeSqr_local = huge(0.0d0)
      directForceCalculations       = 0.0d0
      BHTreeForceCalculations       = 0.0d0

   end subroutine initializeVariables
      
   !------------------------------
   ! Initialize
   !------------------------------   
   subroutine initializeDataArrays
      implicit none

      allocate(body_x   (nrOfBodies)); allocate(body_y       (nrOfBodies))
      allocate(body_vx  (nrOfBodies)); allocate(body_vy      (nrOfBodies))
      allocate(body_ax  (nrOfBodies)); allocate(body_ay      (nrOfBodies))
      allocate(body_jx  (nrOfBodies)); allocate(body_jy      (nrOfBodies))
      allocate(body_mass(nrOfBodies)); allocate(body_mass_inv(nrOfBodies))
      allocate(body_mu  (nrOfBodies))

      ! TODO: arrays for previous accelerations don't need to store all bodies, but only those for the thread
      ! for the multistep integrators we need to remember the accelerations on previous timesteps
      if (integrator .eq. Integrator_Multistep4) then
         allocate(body_axp(3,nrOfBodies))
         allocate(body_ayp(3,nrOfBodies))
      elseif (integrator .eq. Integrator_Multistep6) then
         allocate(body_axp(5,nrOfBodies))
         allocate(body_ayp(5,nrOfBodies))
      elseif (integrator .eq. Integrator_Multistep8) then
         allocate(body_axp(7,nrOfBodies))
         allocate(body_ayp(7,nrOfBodies))
      end if

   end subroutine initializeDataArrays
   
   !------------------------------
   ! Initialize
   !------------------------------   
   subroutine initializeBodies
      
      implicit none
      
      logical                  :: coIsVirtual, objectIsValid, isDefault
      integer                  :: i, j, nrOfInits, nrOfObjects, bodyID, counterStart
      character(len=value_len) :: initString, initType, rotation, formString
      double precision         :: velMax, velMin, velFactorMin, velFactorMax, velX, velY
      double precision         :: centerX, centerY, sizeX, sizeY, radius, radius2
      double precision         :: coMass, minMass, maxMass
      double precision         :: posX, posY, vel, ang, rad, fac1=1.0d0, fac2=1.0d0
      
      nrOfInits = getInputValue("nr_of_inits", 1)
      
      bodyID = 0
      
      ! init the bodies
      do i = 1, nrOfInits
         write (initString, "(A,I0,A)") "init_", i, "_"
         
         initType = getInputValue(trim(initString)//"src", "object")
         
         if (trim(initType) .eq. "object") then
            bodyID = bodyID + 1
            body_x       (bodyID) = getValue(trim(initString)//"x"   , 0.0d0, isDefault)
            body_y       (bodyID) = getValue(trim(initString)//"y"   , 0.0d0, isDefault)
            body_vx      (bodyID) = getValue(trim(initString)//"vx"  , 0.0d0, isDefault)
            body_vy      (bodyID) = getValue(trim(initString)//"vy"  , 0.0d0, isDefault)
            body_mass    (bodyID) = getValue(trim(initString)//"mass", 1.0d0, isDefault)
            body_mass_inv(bodyID) = 1.0d0 / body_mass(bodyID)
            body_mu      (bodyID) = gravitationalConstant * body_mass(bodyID)
         end if
         
         if (trim(initType) .eq. "system") then
            
            nrOfObjects  = getInputValue(trim(initString)//"nr_of_objects"           , 10               )
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
               bodyID                = bodyID + 1
               body_x       (bodyID) = centerX
               body_y       (bodyID) = centerY
               body_vx      (bodyID) = velX
               body_vy      (bodyID) = velY
               body_mass    (bodyID) = coMass
               body_mass_inv(bodyID) = 1.0d0 / coMass
               body_mu      (bodyID) = gravitationalConstant * body_mass(bodyID)
               counterStart = 2      
            end if
            
            do j = counterStart, nrOfObjects
               objectIsValid = .false.

               do while (.not. objectIsValid)
                  posX = getRandomNumber(dble(-radius), dble(radius))
                  posY = getRandomNumber(dble(-radius), dble(radius))
                  
                  if (abs(posX) .lt. 1.0d-14) posX = 1.0d-14
                  
                  rad = sqrt(posX*posX + posY*posY)
                  
                  if ((rad .ge. radius2) .and. (rad .le. radius)) objectIsValid = .true.
               end do

               if (trim(rotation) .eq. "both") then
                  if (getRandomNumber(0.0d0, 2.0d0) .gt. 1.0d0) then
                     fac1 =  1.0d0; fac2 = -1.0d0
                  else
                     fac1 = -1.0d0; fac2 =  1.0d0
                  end if
               end if

               vel = sqrt(gravitationalConstant * coMass / rad) * getRandomNumber(velFactorMin, velFactorMax)
            
               ang = atan(abs(posY/posX))

               if ((posX .lt. 0.0d0) .and. (posY .ge. 0.0d0)) ang = -ang +             PI
               if ((posX .lt. 0.0d0) .and. (posY .lt. 0.0d0)) ang =  ang +             PI
               if ((posX .gt. 0.0d0) .and. (posY .lt. 0.0d0)) ang = -ang + 2.0d0 * PI

               bodyID                = bodyID + 1
               body_x       (bodyID) = centerX + posX
               body_y       (bodyID) = centerY + posY
               body_vx      (bodyID) = velX + fac1 * vel * sin(ang)
               body_vy      (bodyID) = velY + fac2 * vel * cos(ang)
               body_mass    (bodyID) = getRandomNumber(dble(minMass), dble(maxMass))
               body_mass_inv(bodyID) = 1.0d0 / body_mass(bodyID)
               body_mu      (bodyID) = gravitationalConstant * body_mass(bodyID)
            
            end do
         end if
         
         if (trim(initType) .eq. "random") then
         
            nrOfObjects  = getInputValue(trim(initString)//"nr_of_objects"      , 10                    )
            formString   = getValue(trim(initString)//"form"                    , "circle"   , isDefault)
            centerX      = getValue(trim(initString)//"center_x"                , 0.0d0      , isDefault)
            centerY      = getValue(trim(initString)//"center_y"                , 0.0d0      , isDefault)
            velMax       = getValue(trim(initString)//"vel_max"                 , 0.0d0      , isDefault)
            velMin       = getValue(trim(initString)//"vel_min"                 , 0.0d0      , isDefault)
            minMass      = getValue(trim(initString)//"min_mass"                , 1.0d20     , isDefault)
            maxMass      = getValue(trim(initString)//"max_mass"                , 1.0d25     , isDefault)
            
            if (trim(formString) .eq. "circle") then
               radius = getValue(trim(initString)//"radius"                , 1.0d10     , isDefault)
               sizeX  = radius * 2.0d0
               sizeY  = radius * 2.0d0
            elseif (trim(formString) .eq. "rectangle") then
               sizeX  = getValue(trim(initString)//"size_x"                , 1.0d10     , isDefault)
               sizeY  = getValue(trim(initString)//"size_y"                , 1.0d10     , isDefault)
               radius = huge(0.0d0)
            else
               write (*,*) " WARNING: Unkown form parameter for system ", i
               sizeX  = 1d11
               sizeY  = 1d11
               radius = huge(0.0d0)
            end if
            
            do j = 1, nrOfObjects

               ! get position
               objectIsValid = .false.
               do while (.not. objectIsValid)
                  posX = getRandomNumber(-sizeX/2.0d0, sizeX/2.0d0)
                  posY = getRandomNumber(-sizeY/2.0d0, sizeY/2.0d0)

                  rad = sqrt(posX*posX + posY*posY)
                  
                  if (rad .le. radius) objectIsValid = .true.
               end do
               
               ! get velocity
               objectIsValid = .false.
               do while (.not. objectIsValid)
                  velX = getRandomNumber(-velMax, velMax)
                  velY = getRandomNumber(-velMax, velMax)

                  vel = sqrt(velX*velX + velY*velY)
                  
                  if ((vel .ge. velMin) .and. (vel .le. velMax)) objectIsValid = .true.
               end do
               
               bodyID                = bodyID + 1
               body_x       (bodyID) = centerX + posX
               body_y       (bodyID) = centerY + posY
               body_vx      (bodyID) = velX
               body_vy      (bodyID) = velY
               body_mass    (bodyID) = getRandomNumber(dble(minMass), dble(maxMass))
               body_mass_inv(bodyID) = 1.0d0 / body_mass(bodyID)
               body_mu      (bodyID) = gravitationalConstant * body_mass(bodyID)
            end do
         end if
      end do
      
   end subroutine initializeBodies
   
   !------------------------------
   ! Get the number of bodies from the input file
   !------------------------------   
   subroutine getNumberOfBodies
      
      implicit none
      
      logical                  :: isDefault
      integer                  :: i, nrOfInits, nrOfObjects
      character(len=value_len) :: initString, initType

      nrOfInits = getValue("nr_of_inits", 1, isDefault)
      if (isDefault) then
         write (*,*) "Number of inits not specified! STOP"
         call finishIniFileReader()
         !call MPI_FINALIZE (mpiErr) 
         stop
      end if
      
      nrOfBodies = 0
      
      ! init the bodies
      do i = 1, nrOfInits
         write (initString, "(A,I0,A)") "init_", i, "_"
         
         initType = getValue(trim(initString)//"src", "object", isDefault)
         
         if (trim(initType) .eq. "object") then
            nrOfBodies   = nrOfBodies + 1
         end if
         
         if (trim(initType) .eq. "system" .or. trim(initType) .eq. "random") then
            nrOfObjects  = getValue(trim(initString)//"nr_of_objects", 10, isDefault)
            nrOfBodies   = nrOfBodies + nrOfObjects
         end if
      end do
   
   end subroutine getNumberOfBodies
      
   !------------------------------
   ! Initialize
   !------------------------------   
   subroutine selectIntegrator
      implicit none
      
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
            write (*,*) "Unkown integrator selected: ", trim(integratorString), "! Using leapfrog."
            integrator = Integrator_Leapfrog
      end select
   end subroutine selectIntegrator
end module initialize
