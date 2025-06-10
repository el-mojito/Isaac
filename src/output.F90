! ##################################
! OUTPUT ROUTINES
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
   subroutine initializeOutput
      implicit none
      
      call system("mkdir -p "//trim(projectName))
      
      if (outputToSingleFile) then
         open (unit = fid_State    , file = trim(projectName)//"/"//trim(projectName)//".dat"        , status="replace")
      else
         open (unit = fid_Database , file = trim(projectName)//"/"//trim(projectName)//fname_Database, status="replace") 
      end if
      
      if (outputEnergyConservation) then
         open (unit = fid_Energies , file = trim(projectName)//"/"//fname_Energies                   , status="replace")
      end if
      
      if (outputTimesteps) then
         open (unit = fid_Timesteps, file = trim(projectName)//"/"//fname_Timesteps                  , status="replace")
      end if
           
   end subroutine initializeOutput
   
   !------------------------------
   ! Close the output files
   !------------------------------      
   subroutine finishOutput
      implicit none
      
      if (outputToSingleFile) then
         close (fid_State)
      else
         close (fid_Database) 
      end if
      
      if (outputEnergyConservation) then
         close (fid_Energies)
      end if
      
      if (outputTimesteps) then
         close (fid_Timesteps)
      end if
   end subroutine finishOutput
      
   !------------------------------
   ! Write the energy file
   !------------------------------      
   subroutine writeEnergyFile
      implicit none
      
      write (fid_Energies, 1010) t, e_kin, e_pot, e_kin + e_pot, (e_kin + e_pot - e_in) / e_in
      
      1010 format(F18.4,4(1x,E18.10))     !,1x,E18.10,1x,E18.10,1x,E18.10,1x,E18.10)
      
   end subroutine writeEnergyFile
      
   !------------------------------
   ! Write the timestep file
   !------------------------------      
   subroutine writeTimestepFile
      implicit none
      
      write (fid_Timesteps, 1020) t, dt, dt_ratio
      
      1020 format(F18.4,1x,F18.4,1x,F18.4)
   end subroutine writeTimestepFile     
    
   !------------------------------
   ! Write the state file
   !------------------------------
   subroutine writeStateFile
      implicit none
      
      integer            :: i
      character(len=200) :: fname, time
      
      if (.not. outputToSingleFile) then
         write (time , '(F0.2)') t
         
         do while (len_trim(time) .lt. 20)
            time = "0" // trim(time)
         end do
      
         write (fname, '(A,A,A,A)') trim(projectName), "_", trim(time), ".dat"
         
         ! add the entry to the database file, then close it to write the file and re-open for the next write
         ! TODO: unbuffered write?
         write (fid_Database, '(A)'  ) trim(fname)
         close (fid_Database) 
         open  (unit = fid_Database, file = trim(projectName)//"/"//trim(projectName)//fname_Database, &
                status="old", position='append')
         
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


