module integrators

   use parameters
   use bodies
   use time
   use barnes_hut

   implicit none
   
   integer         , private :: i
   double precision, private :: dt_base
   
   contains
   
! ##################################
! INTEGRATOR EULER EXPLICIT
! ##################################
double precision function evolve_euler_explicit()

   ! calculate the accelerations
   evolve_euler_explicit = calculateAccelerations()

   !$omp do
   do i = 1, nrOfBodies
      body_x (i) = body_x (i) + body_vx(i) * dt
      body_y (i) = body_y (i) + body_vy(i) * dt
      
      body_vx(i) = body_vx(i) + body_ax(i) * dt
      body_vy(i) = body_vy(i) + body_ay(i) * dt
   end do
   !$omp end do
   
end function evolve_euler_explicit



! ##################################
! INTEGRATOR LEAPFROG
! ##################################
double precision function evolve_leapfrog()

   !TODO: this prevents a strange compile error...
   i = i
   
   !$omp do
   do i = 1, nrOfBodies
      body_vx(i) = body_vx(i) + 0.5d0 * body_ax(i) * dt
      body_vy(i) = body_vy(i) + 0.5d0 * body_ay(i) * dt
      body_x (i) = body_x (i) +         body_vx(i) * dt
      body_y (i) = body_y (i) +         body_vy(i) * dt
   end do
   !$omp end do

   ! calculate the accelerations
   evolve_leapfrog = calculateAccelerations()

   !$omp do
   do i = 1, nrOfBodies
      body_vx(i) = body_vx(i) + 0.5d0 * body_ax(i) * dt
      body_vy(i) = body_vy(i) + 0.5d0 * body_ay(i) * dt
   end do
   !$omp end do

end function evolve_leapfrog



! ##################################
! INTEGRATOR RUNGE KUTTA 2
! ##################################
double precision function evolve_runge_kutta_2()

end function evolve_runge_kutta_2



! ##################################
! INTEGRATOR RUNGE KUTTA 4
! ##################################
double precision function evolve_runge_kutta_4()

end function evolve_runge_kutta_4



! ##################################
! INTEGRATOR YOSHIDA 4
! ##################################
double precision function evolve_yoshida_4()

   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable

   !$omp barrier
   dt =  1.351207191959657d0 * dt_base
   evolve_yoshida_4 = evolve_leapfrog()
   !$omp barrier
   dt = -1.702414383919315d0 * dt_base
   evolve_yoshida_4 = evolve_leapfrog()
   
   do_dtCalculation = dt_CalculationEnable
   
   !$omp barrier
   dt =  1.351207191959657d0 * dt_base
   evolve_yoshida_4 = evolve_leapfrog()
   
   dt =  dt_base
   
end function evolve_yoshida_4



! ##################################
! INTEGRATOR YOSHIDA 6
! ##################################
double precision function evolve_yoshida_6()

   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable
   
   !$omp barrier 
   dt =  0.784513610477560d0 * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   !$omp barrier
   dt =  0.235573213359357d0 * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   !$omp barrier
   dt = -1.17767998417887d0  * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   
   !$omp barrier
   dt =  1.31518632068391d0  * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   
   !$omp barrier
   dt = -1.17767998417887d0  * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   !$omp barrier
   dt =  0.235573213359357d0 * dt_base
   evolve_yoshida_6 = evolve_leapfrog()
   
   do_dtCalculation = dt_CalculationEnable
   
   !$omp barrier
   dt =  0.784513610477560d0 * dt_base
   evolve_yoshida_6 = evolve_leapfrog()

   dt =  dt_base

end function evolve_yoshida_6



! ##################################
! INTEGRATOR YOSHIDA 8
! ##################################
double precision function evolve_yoshida_8()

   dt_base          = dt
   do_dtCalculation = dt_CalculationDisable

   !$omp barrier
   dt =  0.104242620869991d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt =  0.182020630970714d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt =  0.157739928123617d0     * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt =  0.244002732616735d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt = -0.716989419708120d-2  * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt = -0.244699182370524d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt = -0.161582374150097d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()

   !$omp barrier
   dt = -0.17808286265894516d1 * dt_base
   evolve_yoshida_8 = evolve_leapfrog()

   !$omp barrier
   dt = -0.161582374150097d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()
   !$omp barrier
   dt = -0.244699182370524d1   * dt_base
   evolve_yoshida_8 =  evolve_leapfrog()
   !$omp barrier
   dt = -0.716989419708120d-2  * dt_base
   evolve_yoshida_8 =  evolve_leapfrog()
   !$omp barrier
   dt =  0.244002732616735d1   * dt_base
   evolve_yoshida_8 =  evolve_leapfrog()
   !$omp barrier
   dt =  0.157739928123617d0   * dt_base
   evolve_yoshida_8 =  evolve_leapfrog()
   !$omp barrier
   dt =  0.182020630970714d1   * dt_base
   evolve_yoshida_8 =  evolve_leapfrog()
   
   do_dtCalculation = dt_CalculationEnable
   
   !$omp barrier
   dt =  0.104242620869991d1   * dt_base
   evolve_yoshida_8 = evolve_leapfrog()

   dt =  dt_base

end function evolve_yoshida_8



! ##################################
! REV UP MULTISTEP INTEGRATORS
! ##################################
subroutine rev_up_multistep()
   
   abstract interface
      double precision function evolver ()
      end function evolver
   end interface

   procedure (evolver), pointer  :: evolve => null ()

   integer                       :: i, rev_up_timesteps
   double precision, allocatable :: body_start_x (:), body_start_y (:)
   double precision, allocatable :: body_start_vx(:), body_start_vy(:)
   
   ! TODO: variable timesteps?
       
   allocate (body_start_x (nrOfBodies))
   allocate (body_start_y (nrOfBodies))
   allocate (body_start_vx(nrOfBodies))
   allocate (body_start_vy(nrOfBodies))
   
   ! save starting positions & velocities
   body_start_x  = body_x
   body_start_y  = body_y
   body_start_vx = body_vx
   body_start_vy = body_vy
   
   dt = -dt
   
   if (integrator .eq. Integrator_MultiStep4) then
      rev_up_timesteps =  3
      evolve           => evolve_yoshida_4
   elseif (integrator .eq. Integrator_MultiStep6) then
      rev_up_timesteps =  5
      evolve           => evolve_yoshida_6
   elseif (integrator .eq. Integrator_MultiStep8) then
      rev_up_timesteps =  7
      evolve           => evolve_yoshida_8
   else
      write (*,*) "subroutine rev_up_multistep() !called with unknown multistep integrator! Using MS4."
      integrator       = Integrator_MultiStep4
      rev_up_timesteps =  3
      evolve           => evolve_yoshida_4
   end if
   
   ! rev up
   do i = 1, rev_up_timesteps
      !call evolve()
      body_axp(i,:) = body_ax(1:nrOfBodies)
      body_ayp(i,:) = body_ay(1:nrOfBodies)
   end do
   
   ! reset the system to the starting positions & velocities
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

   double precision :: jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr_local = huge(0.0d0)
   
   ! calculate the accelerations
   !call calculateAccelerations()
   
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
   
   body_axp(1,:) = body_ax(1:nrOfBodies)
   body_ayp(1,:) = body_ay(1:nrOfBodies)

end subroutine evolve_multistep_4



! ##################################
! INTEGRATOR MULTISTEP 4
! ##################################
subroutine evolve_multistep_6()

   !double precision :: jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr_local = huge(0.0d0)
   
end subroutine evolve_multistep_6



! ##################################
! INTEGRATOR MULTISTEP 8
! ##################################
subroutine evolve_multistep_8()

   !double precision :: jdt(2), sdt2(2), cdt3(2)
      
   ! TODO: variable timesteps?   
   
   minimumCollisionTimeSqr_local = huge(0.0d0)

end subroutine evolve_multistep_8

end module integrators
