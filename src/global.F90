!##############################
! Parameter
!##############################
module parameters
   
   implicit none
   
   integer         , parameter :: Integrator_EulerExplicit      = 1
   integer         , parameter :: Integrator_Leapfrog           = 2
   integer         , parameter :: Integrator_RungeKutta2        = 3
   integer         , parameter :: Integrator_RungeKutta4        = 4
   integer         , parameter :: Integrator_RungeKutta6        = 5
   integer         , parameter :: Integrator_Yoshida4           = 6
   integer         , parameter :: Integrator_Yoshida6           = 7
   integer         , parameter :: Integrator_Yoshida8           = 8
   integer         , parameter :: Integrator_MultiStep4         = 9
   integer         , parameter :: Integrator_MultiStep6         = 10
   integer         , parameter :: Integrator_MultiStep8         = 11
   integer         , parameter :: Integrator_PredictorCorrector = 12
   integer         , parameter :: Integrator_Hermite4           = 13
   integer         , parameter :: Integrator_Hermite6           = 14
   integer         , parameter :: Integrator_Hermite8           = 15
   
   integer         , parameter :: Force_Direct                  = 1
   integer         , parameter :: Force_BarnesHut               = 2
   
   integer         , parameter :: Timestep_Constant             = 1
   integer         , parameter :: Timestep_Variable             = 2
      
   double precision, parameter :: PI                            = 2.0d0*acos(0.0d0)
end module parameters


!##############################
! Global
!##############################
module global
   use ini_file_reader
   
   implicit none
      
   logical                  :: doRemoveRunaways, passParameterChange
   integer                  :: integrator, forceCalculation
   integer                  :: randomSeeds, randomSeed(1)
   integer                  :: removeRunawaysIterations
   double precision         :: gravitationalConstant
   double precision         :: removeRunawaysDistanceFromCom
   double precision         :: maxDistance
   double precision         :: directForceCalculations, BHTreeForceCalculations
   character(len=value_len) :: projectName
   character(len=value_len) :: integratorString, forceString
   character(len=100)       :: inputFile
   
   contains
         
   !------------------------------
   ! Get a random number between two limits
   !------------------------------
   double precision function getRandomNumber(lower, upper)
      implicit none
   
      double precision, intent(in) :: lower, upper
      real                         :: rnd
   
      call random_number(rnd)
   
      getRandomNumber = lower + dble(rnd) * (upper - lower)
   end function getRandomNumber

end module global


!##############################
! Bodies
!##############################
module bodies
   implicit none
   save

   integer                       :: nrOfBodies, bodiesRemoved
   double precision, allocatable :: body_x (:)   , body_y (:)                  ! body coordinates
   double precision, allocatable :: body_vx(:)   , body_vy(:)                  ! body velocities
   double precision, allocatable :: body_ax(:)   , body_ay(:)                  ! body accelerations
   double precision, allocatable :: body_axp(:,:), body_ayp(:,:)               ! body accelerations on previous timesteps for multistep integrator
   double precision, allocatable :: body_jx(:)   , body_jy(:)                  ! body jerks
   double precision, allocatable :: body_mass(:) , body_mass_inv(:)            ! body mass, inverted mass
   double precision, allocatable :: body_mu(:)                                 ! gravitational parameter, mu = gravitational constant * body mass

   contains
     
   !------------------------------
   ! Calculate the maximum distance one of the bodies has on any axis from the center
   !------------------------------
   double precision function getMaxXYDistanceFromCenter()
      implicit none
      
      integer          :: i
      double precision :: maxDistance
      
      maxDistance = 0.0d0
      
      !TODO: improve scheduling
      !$omp parallel do reduction(max : maxDistance) schedule(dynamic)
      do i = 1, nrOfBodies
         maxDistance = max(maxDistance, abs(body_x(i)), abs(body_y(i)))
      
         !if (abs(body_x(i)) .gt. maxDistance) maxDistance = abs(body_x(i))
         !if (abs(body_y(i)) .gt. maxDistance) maxDistance = abs(body_y(i))
      end do
      !$omp end parallel do
      
      getMaxXYDistanceFromCenter = maxDistance
      
   end function getMaxXYDistanceFromCenter
end module bodies


!##############################
! MPI PARALLEL
!##############################
!module parallel_mpi
!   use mpi
!   use bodies

!   implicit none
!   
!   integer, parameter   :: mpiMaster = 0
!   
!   integer              :: mpiThreads, mpiID, mpiErr
!   integer, allocatable :: mpiCounterStartList(:), mpiCounterEndList(:), mpiBodyCount(:), mpiOffset(:)

!   contains
!   
!   !------------------------------
!   ! Synchronize the bodies among all threads
!   !------------------------------
!   subroutine mpiSynchronizeBodies()
!      call MPI_ALLGATHERV(MPI_IN_PLACE, nrOfBodies             , MPI_DOUBLE_PRECISION, &
!                          body_x      , mpiBodyCount, mpiOffset, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
!      call MPI_ALLGATHERV(MPI_IN_PLACE, nrOfBodies             , MPI_DOUBLE_PRECISION, &
!                          body_y      , mpiBodyCount, mpiOffset, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
!      call MPI_ALLGATHERV(MPI_IN_PLACE, nrOfBodies             , MPI_DOUBLE_PRECISION, &
!                          body_vx     , mpiBodyCount, mpiOffset, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
!      call MPI_ALLGATHERV(MPI_IN_PLACE, nrOfBodies             , MPI_DOUBLE_PRECISION, &
!                          body_vy     , mpiBodyCount, mpiOffset, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
!   end subroutine
!   
!   !------------------------------
!   ! Allocate the necessary arrays to divide the bodies among all threads
!   !------------------------------
!   subroutine mpiInitializeArrays()
!      allocate(mpiCounterStartList(mpiThreads))
!      allocate(mpiCounterEndList  (mpiThreads))
!      allocate(mpiBodyCount       (mpiThreads))
!      allocate(mpiOffset          (mpiThreads))
!   end subroutine
!end module parallel_mpi


!##############################
! OMP PARALLEL
!##############################
module parallel_omp
   implicit none
   
   integer              :: ompThreads

end module parallel_omp


!##############################
! Time
!##############################
module time
   use ini_file_reader
   
   implicit none
      
   logical                  :: reachedEndTime, do_dtCalculation, dt_CalculationEnable, dt_CalculationDisable
   integer(kind=16)         :: iteration
   integer                  :: timestepMode
   real(kind=4)             :: timearray(2)
   double precision         :: startTimeInit, startTimeCalc, startTimeUser
   double precision         ::   endTimeInit,   endTimeCalc, timeNow
   double precision         :: t, t_end, dt, dt_min, dt_max, dt_old, dt_ratio, dt_user
   double precision         :: dt_factor
   double precision         :: minimumCollisionTimeSqr, minimumCollisionTimeSqr_local
   character(len=value_len) :: timestepModeString
   
   contains
   
   !------------------------------
   ! This subroutine calculates the time step from the minimum collision time
   ! and keeps track of the maximum and minimum time step used
   !------------------------------   
   !subroutine calculateTimeStep
   !
   !   dt_old   = dt  
   !   if (timestepMode .eq. Timestep_Variable) then
   !      dt = sqrt(minimumCollisionTimeSqr) * dt_factor
   !   end if
   !      
   !   if (outputExactTimes) then
   !        if (t+dt > t_output) dt = t_output - t
   !        if (t+dt > t_end)    dt = t_end    - t
   !   end if
   !   dt_ratio = dt / dt_old
   !      
   !   if (dt .gt. dt_max) dt_max = dt
   !   if (dt .lt. dt_min) dt_min = dt
   !   
   !end subroutine calculateTimeStep
end module time


!##############################
! Output
!##############################
module output_data
   
   implicit none
  
   integer          , parameter :: fid_Input     = 50
   integer          , parameter :: fid_Energies  = 51
   integer          , parameter :: fid_Database  = 52
   integer          , parameter :: fid_Timesteps = 53
   integer          , parameter :: fid_State     = 54
   
   character(len=12), parameter :: fname_Energies  = 'energies.dat'
   character(len=15), parameter :: fname_Database  = '_database.visit'
   character(len=13), parameter :: fname_Timesteps = 'timesteps.dat'
   
   !logical            :: outputEnabled
   logical            :: outputToSingleFile, outputEnergyConservation, outputTimesteps, outputSilentMode, outputExactTimes
   integer            :: percentage, lastPercentageWritten
   double precision   :: t_output, dt_output 
end module output_data


!##############################
! Energy
!##############################
module energy
   !use parallel_mpi
   
   implicit none
   
   double precision ::     e_kin,     e_pot,     e_in,     e_out
   
   contains
   
   subroutine calculateEnergies
      use global
      use bodies
      
      implicit none
      
      integer          :: i, j
      
      e_kin = 0.0d0
      e_pot = 0.0d0
      
      !TODO: improve scheduling
      !$omp parallel do private(j) reduction(+ : e_kin, e_pot) schedule(dynamic)
      do i = 1, nrOfBodies
         e_kin = e_kin + 0.5d0 * body_mass(i) * (body_vx(i)*body_vx(i) + body_vy(i)*body_vy(i))
         do j = i+1, nrOfBodies
            e_pot = e_pot - gravitationalConstant * (body_mass(i)*body_mass(j)) / & 
                         sqrt((body_x(i)-body_x(j))**2 + (body_y(i)-body_y(j))**2)
         end do
      end do
      !$omp end parallel do
      
   end subroutine calculateEnergies
end module energy


!##############################
! Quads
!##############################
module quads
   
   implicit none

   integer         , parameter :: QERROR=0, QNW=1, QNE=2, QSW=3, QSE=4
   double precision, parameter :: dir_north = 1.0d0, dir_south = -1.0d0
   double precision, parameter :: dir_east  = 1.0d0, dir_west  = -1.0d0

   type quad
      double precision :: center_x, center_y
      double precision :: length, lengthSqr
   end type quad
   
   contains
   
   !------------------------------
   ! Check if a body lays within a quad
   !------------------------------
   pure logical function containsBody(q, bid)
      use bodies
      implicit none
      
      type(quad), intent(in) :: q
      integer,    intent(in) :: bid
      double precision       :: lh
      
      lh = q%length / 2.0d0
      
      if ( q%center_x - lh <= body_x(bid) .and. &
           q%center_x + lh >= body_x(bid) .and. &
           q%center_y - lh <= body_y(bid) .and. &
           q%center_y + lh >= body_y(bid) ) then
         containsBody=.true.
      else
         containsBody=.false.
      end if
   end function containsBody
      
   !------------------------------
   ! Get an internal quad
   !------------------------------
   pure type(quad) function getInternalQuad(q, dir_x, dir_y)
      implicit none
      type(quad),       intent(in) :: q
      double precision, intent(in) :: dir_x, dir_y
      double precision             :: len4
      
      len4 = q%length / 4.0d0
      getInternalQuad%center_x  = q%center_x + dir_x * len4
      getInternalQuad%center_y  = q%center_y + dir_y * len4
      getInternalQuad%length    = len4 * 2.0d0
      getInternalQuad%lengthSqr = getInternalQuad%length * getInternalQuad%length
   end function getInternalQuad
   
      
   !------------------------------
   ! Get the sub-quad in which the body is
   !------------------------------
   integer function getSubQuadrant(q, bid)
      use bodies
      use time
      implicit none
      
      type(quad), intent(in) :: q
      integer,    intent(in) :: bid
      integer                :: res
      
      res = QERROR
      
      if ((body_x(bid) .le. q%center_x) .and. (body_y(bid) .ge. q%center_y)) res = QNW
      if ((body_x(bid) .gt. q%center_x) .and. (body_y(bid) .ge. q%center_y)) res = QNE
      if ((body_x(bid) .le. q%center_x) .and. (body_y(bid) .lt. q%center_y)) res = QSW
      if ((body_x(bid) .gt. q%center_x) .and. (body_y(bid) .lt. q%center_y)) res = QSE

      getSubQuadrant = res
   end function getSubQuadrant
end module quads


!##############################
! Barnes-Hut
!##############################
module barnes_hut
   !use mpi
   !use parallel_mpi
   use omp_lib
   use global
   use bodies
   use quads
   use time
   
   implicit none
   
   integer          :: nodes
   integer          :: idCounter
   double precision :: theta, thetaSqr

   type barnes_hut_tree_node
      integer          :: id, body_id
      integer          :: nrOfBodies
      double precision :: com_x, com_y
      double precision :: mass
      type(quad)       :: quad
      type(barnes_hut_tree_node), pointer :: NW=>null(), NE=>null(), SW=>null(), SE=>null()
   end type barnes_hut_tree_node
   
   type (barnes_hut_tree_node) :: BH_Tree
   
   contains
   
   !------------------------------
   ! Build the BH Tree
   !------------------------------
   subroutine buildBHTree()
      implicit none
      
      integer :: i

      call resetBHTree()
      
      maxDistance = getMaxXYDistanceFromCenter()
      call initBHTreeQuad(2.01d0 * maxDistance)
      
      do i = 1, nrOfBodies
        call insertBodyInBHTree(BH_Tree, i)
      end do
         
      call calculateBHTreeMassDistribution()
     
   end subroutine buildBHTree
      
   !------------------------------
   ! Reset the BH Tree
   !------------------------------
   subroutine resetBHTree()
      implicit none

      nodes                   = 0
      idCounter               = 0
      BH_Tree%id              = 0
      BH_Tree%body_id         = 0
      BH_Tree%nrOfBodies      = 0
      BH_Tree%mass            = 0.0d0
      BH_Tree%com_x           = 0.0d0
      BH_Tree%com_y           = 0.0d0
      
      if (associated(BH_Tree%NW)) then
         call resetBHTreeNode(BH_Tree%NW); deallocate(BH_Tree%NW); nullify(BH_Tree%NW)
      end if
      if (associated(BH_Tree%NE)) then
         call resetBHTreeNode(BH_Tree%NE); deallocate(BH_Tree%NE); nullify(BH_Tree%NE)
      end if
      if (associated(BH_Tree%SW)) then
         call resetBHTreeNode(BH_Tree%SW); deallocate(BH_Tree%SW); nullify(BH_Tree%SW)
      end if
      if (associated(BH_Tree%SE)) then
         call resetBHTreeNode(BH_Tree%SE); deallocate(BH_Tree%SE); nullify(BH_Tree%SE)
      end if
     
   end subroutine resetBHTree
   
   !------------------------------
   ! Reset a node
   !------------------------------   
   recursive subroutine resetBHTreeNode(bht)
      implicit none

      type(barnes_hut_tree_node), target :: bht
      
      if (associated(bht%NW)) then
         call resetBHTreeNode(bht%NW); deallocate(bht%NW); nullify(bht%NW)
      end if
      if (associated(bht%NE)) then
         call resetBHTreeNode(bht%NE); deallocate(bht%NE); nullify(bht%NE)
      end if
      if (associated(bht%SW)) then
         call resetBHTreeNode(bht%SW); deallocate(bht%SW); nullify(bht%SW)
      end if
      if (associated(bht%SE)) then
         call resetBHTreeNode(bht%SE); deallocate(bht%SE); nullify(bht%SE)
      end if

   end subroutine resetBHTreeNode   
   
   !------------------------------
   ! Initalize the main quad
   !------------------------------   
   subroutine initBHTreeQuad(length)
      implicit none
      
      double precision, intent(in) :: length
      
      BH_Tree%quad%center_x  = 0.0d0
      BH_Tree%quad%center_y  = 0.0d0
      BH_Tree%quad%length    = length
      BH_Tree%quad%lengthSqr = length * length

   end subroutine initBHTreeQuad      
   
   !------------------------------
   ! Initalize a BH Tree node
   !------------------------------   
   subroutine initBHTreeNode(bht, q)
      implicit none
      
      type(quad)                , intent(in)    :: q
      type(barnes_hut_tree_node), intent(inout) :: bht
      
      idCounter = idCounter + 1
      
      bht%id         = idCounter
      bht%body_id    = 0
      bht%nrOfBodies = 0
      bht%quad       = q
      bht%mass       = 0.0d0
      bht%com_x      = 0.0d0
      bht%com_y      = 0.0d0
      nullify(bht%NW)
      nullify(bht%NE) 
      nullify(bht%SW)
      nullify(bht%SE)
    
   end subroutine initBHTreeNode   
   
   !------------------------------
   ! Insert a body in the BH Tree
   !------------------------------   
   subroutine insertBodyInBHTree(bht, bid)
      implicit none

      type(barnes_hut_tree_node), intent(inout) :: bht      ! BH Treee node
      integer                   , intent(in)    :: bid      ! body id
      
      if (bht%nrOfBodies .gt. 1) then
         call moveBodyToSubNode(bht, bid)
      elseif (bht%nrOfBodies .eq. 1) then
         call moveBodyToSubNode(bht, bht%body_id)
         call moveBodyToSubNode(bht, bid)   
      else
         bht%body_id = bid
      end if
      
      bht%nrOfBodies = bht%nrOfBodies + 1

   end subroutine insertBodyInBHTree   
       
   !------------------------------
   ! Move a body one node deeper
   !------------------------------   
   subroutine moveBodyToSubNode(bht, bid)
      implicit none

      type(barnes_hut_tree_node), intent(inout) :: bht      ! BH Treee node
      integer                   , intent(in)    :: bid      ! body id
      type(barnes_hut_tree_node), pointer       :: subnode
      integer                                   :: quadrant
      
      quadrant = getSubQuadrant(bht%quad, bid)
      
      select case (quadrant)
         case (QNW)
            if (.not. associated(bht%NW)) then
               allocate(bht%NW)
               call initBHTreeNode(bht%NW, getInternalQuad(bht%quad, dir_west, dir_north))
               nodes = nodes + 1
            end if
            subnode => bht%NW
         case (QNE) 
            if (.not. associated(bht%NE)) then
               allocate(bht%NE)
              call initBHTreeNode(bht%NE, getInternalQuad(bht%quad, dir_east, dir_north))
               nodes = nodes + 1
            end if
            subnode => bht%NE
         case (QSW) 
            if (.not. associated(bht%SW)) then
               allocate(bht%SW)
               call initBHTreeNode(bht%SW, getInternalQuad(bht%quad, dir_west, dir_south))
               nodes = nodes + 1
            end if
            subnode => bht%SW
         case (QSE)
            if (.not. associated(bht%SE)) then
               allocate(bht%SE)
               call initBHTreeNode(bht%SE, getInternalQuad(bht%quad, dir_east, dir_south))
               nodes = nodes + 1
            end if
            subnode => bht%SE
         case default
            write (*,*) "Error in moveBodyToSubNode(): could not find quad for body ", bid
            subnode => bht%NW
      end select
      
      call insertBodyInBHTree(subnode, bid) 

   end subroutine moveBodyToSubNode 
         
   !------------------------------
   ! Calculate the BH Tree mass distribution
   !------------------------------   
   subroutine calculateBHTreeMassDistribution()
      implicit none
      
      call calculateNodeMassDistribution(BH_Tree)
      
   end subroutine calculateBHTreeMassDistribution
   
   !------------------------------
   ! Calculate the mass distribution of a node
   !------------------------------   
   recursive subroutine calculateNodeMassDistribution(bht)
      implicit none
      
      type(barnes_hut_tree_node), intent(inout) :: bht
      double precision                          :: mass_inv
      
      if (bht%nrOfBodies .eq. 1) then
         bht%mass  = body_mass(bht%body_id)
         bht%com_x = body_x   (bht%body_id)
         bht%com_y = body_y   (bht%body_id)
      else
         if (associated(bht%NW)) then
            call calculateNodeMassDistribution(bht%NW)
            bht%mass  = bht%mass  + bht%NW%mass
            bht%com_x = bht%com_x + bht%NW%mass * bht%NW%com_x
            bht%com_y = bht%com_y + bht%NW%mass * bht%NW%com_y
         end if
         if (associated(bht%NE)) then
            call calculateNodeMassDistribution(bht%NE)
            bht%mass  = bht%mass  + bht%NE%mass
            bht%com_x = bht%com_x + bht%NE%mass * bht%NE%com_x
            bht%com_y = bht%com_y + bht%NE%mass * bht%NE%com_y
         end if
         if (associated(bht%SW)) then
            call calculateNodeMassDistribution(bht%SW)
            bht%mass  = bht%mass  + bht%SW%mass
            bht%com_x = bht%com_x + bht%SW%mass * bht%SW%com_x
            bht%com_y = bht%com_y + bht%SW%mass * bht%SW%com_y
         end if
         if (associated(bht%SE)) then
            call calculateNodeMassDistribution(bht%SE)
            bht%mass  = bht%mass  + bht%SE%mass
            bht%com_x = bht%com_x + bht%SE%mass * bht%SE%com_x
            bht%com_y = bht%com_y + bht%SE%mass * bht%SE%com_y
         end if 
         
         mass_inv = 1.0d0 / bht%mass
         bht%com_x = bht%com_x * mass_inv
         bht%com_y = bht%com_y * mass_inv
         
         !bht%com_x = bht%com_x / bht%mass
         !bht%com_y = bht%com_y / bht%mass
         
      end if
      
   end subroutine calculateNodeMassDistribution
      
   !------------------------------
   ! Calculate all accelerations
   !------------------------------
   double precision function calculateAccelerations()
      implicit none
      
      integer                       :: i
      double precision              :: f(2)
      double precision              :: colTimeSqr, minColTimeSqr_local
      
      ! if we calculate the accelerations the system has changed (or else we would not need to calculate them)
      ! that means we have to synchronize first and then build a new BH Tree
      ! call mpiSynchronizeBodies()
      
      !$omp single
      call buildBHTree()
      !$omp end single

      minColTimeSqr_local = huge(0.0d0)
      
      !TODO: improve scheduling
      !$omp do schedule(dynamic)
      do i = 1, nrOfBodies
         colTimeSqr = huge(0.0d0)
         call calculateForce(BH_Tree, i, f, colTimeSqr)
         
         minColTimeSqr_local = min(minColTimeSqr_local, colTimeSqr)
         
         !if (colTimeSqr .lt. minColTimeSqr_local) then
         !   minColTimeSqr_local = colTimeSqr
         !   write (*,*) "new minColTime for bid=", i
         !   write (*,*) "            minColTime=", sqrt(colTimeSqr)
         !end if
          
         body_ax(i) = f(1) * body_mass_inv(i)
         body_ay(i) = f(2) * body_mass_inv(i)
      end do
      !$omp end do

      calculateAccelerations = minColTimeSqr_local
      
   end function calculateAccelerations
   
   !------------------------------
   ! Calculate all forces on a body from a tree node
   !------------------------------      
   recursive subroutine calculateForce(bht, bid, f, colTimeSqr)
      implicit none
      
      integer                   , intent(in)    :: bid
      type(barnes_hut_tree_node), intent(in)    :: bht 
      double precision          , intent(out)   :: f (2)
      double precision          , intent(inout) :: colTimeSqr
      double precision                          :: fr(2)
      double precision                          :: r (2), rsqr
      
      f = 0.0d0
       
      if (bht%nrOfBodies .eq. 1) then
         if (bht%body_id .ne. bid) then
            call calculateForceBetweenObjects(bid, body_x(bht%body_id), body_y(bht%body_id), &
                                                   body_mass(bht%body_id), .true., f, &
                                                   bid2=bht%body_id, colTimeSqr=colTimeSqr)
         end if
      else
         r(1) = bht%com_x - body_x(bid)
         r(2) = bht%com_y - body_y(bid)
         
         rsqr = r(1)*r(1) + r(2)*r(2)         ! TODO: idea: calculate before the if? has to be calculated anyway

         if ((bht%quad%lengthSqr / rsqr) < thetaSqr) then
            call calculateForceBetweenObjects(bid, bht%com_x, bht%com_y, bht%mass, .false., &
                                              f, r=r, rsqr=rsqr)
         else
            if (associated(bht%NW)) then
               call calculateForce(bht%NW, bid, fr, colTimeSqr); f = f + fr
            end if
            if (associated(bht%NE)) then
               call calculateForce(bht%NE, bid, fr, colTimeSqr); f = f + fr
            end if
            if (associated(bht%SW)) then
               call calculateForce(bht%SW, bid, fr, colTimeSqr); f = f + fr
            end if
            if (associated(bht%SE)) then
               call calculateForce(bht%SE, bid, fr, colTimeSqr); f = f + fr
            end if
         end if
      end if
      
   end subroutine calculateForce

   !------------------------------
   ! Calculate the force from the object data "bx, by, bz, bm" on the body with id "bid"
   ! bid        = body id to calculate the force on
   ! bx, by     = coordinates of the body which exerts the force
   ! bm         = mass of the body which exerts the force
   ! d          = direct calculation between to bodies (.true.) or between a body and a node (.false.)?
   ! f          = calculated force
   ! r, rsqr    = (optional) distance x, y and distance squared, this was already calculated for body <-> node
   ! bid2       = (optional) body id of the second body, if it is a direct calculation (for the time step calculation)
   ! colTimeSqr = (optional) square collision time for direct force calculation (distance^2 / speed difference^2)
   !------------------------------      
   subroutine calculateForceBetweenObjects(bid, bx, by, bm, d, f, r, rsqr, bid2, colTimeSqr)
      implicit none
      logical         , intent(in)              :: d
      integer         , intent(in)              :: bid
      double precision, intent(in)              :: bx, by, bm
      integer         , intent(in) , optional   :: bid2
      double precision, intent(in) , optional   :: r(2), rsqr
      double precision, intent(out)             :: f(2)
      double precision, intent(inout), optional :: colTimeSqr
      double precision                          :: r21(2), r3
      double precision                          :: v(2)
      double precision                          :: f1

      if (.not. present(r) .or. .not. present(rsqr)) then
         r21(1) = bx - body_x(bid)
         r21(2) = by - body_y(bid)
         r3     = r21(1)*r21(1) + r21(2)*r21(2)                      ! only r² so far
      else
         r21    = r
         r3     = rsqr                                               ! only r² so far
      end if

      !f1 = gravitationalConstant * body_mass(bid) * bm
      f1 = body_mu(bid) * bm                                        ! mu = gravitational constant * body mass
      
      if (d) then
         if (do_dtCalculation) then
            v(1) = body_vx(bid) - body_vx(bid2)
            v(2) = body_vy(bid) - body_vy(bid2)
            
            colTimeSqr = min(colTimeSqr, r3 / (v(1)*v(1) + v(2)*v(2)))

         end if

         !directForceCalculations = directForceCalculations + 1.0d0
      else
         !BHTreeForceCalculations = BHTreeForceCalculations + 1.0d0
      end if
      
      r3 = r3 * sqrt(r3)
      
      f1 = f1 / r3
      
      f = f1 * r21
   end subroutine calculateForceBetweenObjects

end module barnes_hut
