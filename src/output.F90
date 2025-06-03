! ##################################
! MAIN LOOP BARNES HUT
! ##################################
module output
   use global
   use bodies
   use output_data
   use energy
   use time
   
   implicit none
   
   contains
   
   !------------------------------
   ! Initalize the output files
   !------------------------------   
   subroutine initOutput
      implicit none
      
      call system("mkdir -p "//trim(projectName))
      
      open (unit = fid_Energies , file = trim(projectName)//"/"//fname_Energies                    , status="replace")
      open (unit = fid_Timesteps, file = trim(projectName)//"/"//fname_Timesteps                   , status="replace")
      open (unit = fid_Database , file = trim(projectName)//"/"//trim(projectName)//fname_Database , status="replace") 
      
      if (outputToSingleFile) open (unit = fid_State, file = trim(projectName)//"/"//trim(projectName)//".dat", status="replace")
      
   end subroutine initOutput
   
   !------------------------------
   ! Close the output files
   !------------------------------      
   subroutine finishOutput
      implicit none
      
      close (fid_Energies)
      close (fid_Timesteps)
      close (fid_Database) 
      if (outputToSingleFile) close (fid_State)
      
   end subroutine finishOutput
      
   !------------------------------
   ! Write the energy file
   !------------------------------      
   subroutine writeEnergyFile
      implicit none
      
      write (fid_Energies, 1010) t, e_kin, e_pot, e_kin + e_pot, (e_kin + e_pot - e_in) / e_in
      
      1010 format(F15.4,4(1x,E18.10))     !,1x,E18.10,1x,E18.10,1x,E18.10,1x,E18.10)
      
   end subroutine writeEnergyFile
      
   !------------------------------
   ! Write the timestep file
   !------------------------------      
   subroutine writeTimestepFile
      implicit none
            
      write (fid_Timesteps, 1020) dt, dt_ratio
      
      1020 format(F13.4,1x,F13.4)
   end subroutine writeTimestepFile     
    
   !------------------------------
   ! Write the state file
   !------------------------------
   subroutine writeStateFile
      implicit none
      
      integer            :: i
      character(len=200) :: fname
      
      if (.not. outputToSingleFile) then      
         write (fname, '(A,A,F0.2,A)') trim(projectName), "_", t, ".dat"
         
         ! add the entry to the database file
         write (fid_Database, '(A)'  ) trim(fname)
         
         ! write the state file
         open (unit = fid_State, file = trim(projectName)//"/"//trim(fname), status="replace")
         
         write (fid_State, '(A,A,A)' ) "TITLE = """, trim(projectName), """"
         write (fid_State, '(A)'     ) "VARIABLES = ""x[m]"" ,""y[m]"" ,""v_x[m/s]"" ,""v_y[m/s]"" ,""mass[kg]"" ,""id"""
         write (fid_State, '(A)'     ) "ZONE T=""BodyData"""
         write (fid_State, '(A,I0,A)') "I=", nrOfBodies, ", J=1, K=1, F=POINT"
      end if
      
      do i = 1, nrOfBodies
         ! add "buffered='yes'" for ifort
         write (fid_State, *) body_x(i), body_y(i), body_vx(i), body_vy(i), body_mass(i), i
      end do
      
      if (.not. outputToSingleFile) close (fid_State)
   end subroutine writeStateFile

end module output


