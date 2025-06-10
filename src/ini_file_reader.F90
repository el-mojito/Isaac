module ini_file_reader
   use strings
   
   implicit none
   
   integer, parameter                             :: param_len       = 150
   integer, parameter                             :: value_len       = 250

   character(len=param_len), allocatable, private :: ini_param(:)
   character(len=value_len), allocatable, private :: ini_value(:)
   
   logical                              , private :: is_initialized = .false.
   integer                              , private :: file_id
   integer                              , private :: entries_found = 0, current_entry = 0, max_entries = 0
   
   
   interface getValue        
      module procedure getValue_sp     ! single precision
      module procedure getValue_dp     ! double precision
      module procedure getValue_i      ! integer
      module procedure getValue_st     ! string
      module procedure getValue_lo     ! logical
   end interface

   contains
   
   !****************************
   ! 
   !****************************
   subroutine initIniFileReader(fid, mx)
      implicit none
      
      integer, intent(in) :: fid, mx
      
      file_id       = fid
      entries_found = 0
      current_entry = 0
      max_entries   = mx
      
      allocate (ini_param(max_entries))
      allocate (ini_value(max_entries))
      
      is_initialized = .true.
      
   end subroutine initIniFileReader
   
   !****************************
   ! 
   !****************************   
   subroutine finishIniFileReader
      implicit none
      
      deallocate (ini_param)
      deallocate (ini_value)
      
      is_initialized = .false.
      
   end subroutine finishIniFileReader
      
   !****************************
   ! 
   !****************************  
   subroutine analyzeIniFile(ret, entries)
      implicit none
      
      integer, intent(out)                 :: ret
      integer, intent(out), optional       :: entries
      integer                              :: ios, nargs
      character(len=param_len+value_len+3) :: str                 ! +3 to account for possible " = "
      character(len=param_len+value_len+3) :: args(2)             ! +3 to account for possible " = "
      
      ret = 0
      if (present(entries)) entries = 0
            
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Error while analyzing ini file. Reader not initialized!"
         ret = -1
         return
      end if
            
      ios = 0
      do while (ios .ge. 0)
         if (current_entry .eq. max_entries) then
            write (*, "(2A)") "IniFileReader: Error while analyzing ini file. Too many entries in the file.", &
                              "Consider calling initIniFileReader() with a larger number for max_entries!"
            ret = -1
            return
         end if
      
         call readline(file_id, str, ios)
         
         if (ios .lt. 0) return
         
         if (ios .gt. 0) then
            write (*, "(A,I4,2A)") "IniFileReader: Error while analyzing ini file in call to ""readline()"", ios=", &
                                   ios, ". Line: ", trim(str)
            ret = -1
            return
         end if
         
         call parse(str, "=", args, nargs)
         
         if (nargs .gt. 2) then
            write (*, "(3A)") "IniFileReader: Error while analyzing ini file. More than one ""="" found in one line: """, &
                             trim(str), """"
            ret = -1
            return
         end if
         
         if (nargs .eq. 1) then
            args(2) = ""
         end if
         
         entries_found = entries_found + 1
         current_entry = current_entry + 1
         if (present(entries)) entries = entries_found
         
         ini_param(current_entry) = trim(args(1))
         ini_value(current_entry) = trim(args(2))
         
      end do
     
   end subroutine analyzeIniFile
      
   !****************************
   ! 
   !****************************     
   subroutine printIniFile()
      implicit none
      
      integer :: i
      
      do i = 1, entries_found
         write (*,*) "param: ", trim(ini_param(i))
         write (*,*) "value: ", trim(ini_value(i))
      end do
   end subroutine printIniFile
      
   !****************************
   ! 
   !****************************  
   real function getValue_sp(param, def, isdef)
      implicit none
      
      character(len=*), intent(in)  :: param
      real            , intent(in)  :: def
      logical         , intent(out) :: isdef
      integer                       :: i, ios
      real                          :: val
      
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Cannot get value. Reader not initialized! Assigning default value."
         getValue_sp = def
         return
      end if
      
      getValue_sp = def
      isdef       = .true.
      
      do i = 1, entries_found
         if (trim(ini_param(i)) .eq. trim(param)) then
            call value(ini_value(i), val, ios)
            if (ios .eq. 0) then
               getValue_sp = val
               isdef       = .false.
               return
            else
               getValue_sp = def
               return
            end if
         end if
      end do
   end function getValue_sp
      
   !****************************
   ! 
   !****************************  
   double precision function getValue_dp(param, def, isdef)
      implicit none
      
      character(len=*), intent(in)  :: param 
      double precision, intent(in)  :: def
      logical         , intent(out) :: isdef
      integer                       :: i, ios
      double precision              :: val
      
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Cannot get value. Reader not initialized! Assigning default value."
         getValue_dp = def
         return
      end if
      
      getValue_dp = def
      isdef       = .true.
      
      do i = 1, entries_found
         if (trim(ini_param(i)) .eq. trim(param)) then
            call value(ini_value(i), val, ios)
            if (ios .eq. 0) then
               getValue_dp = val
               isdef       = .false.
               return
            else
               getValue_dp = def
               return
            end if
         end if
      end do
   end function getValue_dp
           
   !****************************
   ! 
   !****************************  
   integer function getValue_i(param, def, isdef)
      implicit none
      
      character(len=*), intent(in)  :: param
      integer         , intent(in)  :: def
      logical         , intent(out) :: isdef
      integer                       :: i, ios, val
      
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Cannot get value. Reader not initialized! Assigning default value."
         getValue_i = def
         return
      end if
      
      getValue_i = def
      isdef      = .true.
      
      do i = 1, entries_found
         if (trim(ini_param(i)) .eq. trim(param)) then
            call value(ini_value(i), val, ios)
            if (ios .eq. 0) then
               getValue_i = val
               isdef      = .false.
               return
            else
               getValue_i = def
               return
            end if
         end if
      end do
   end function getValue_i
         
   !****************************
   ! 
   !****************************  
   character(len=value_len) function getValue_st(param, def, isdef)
      implicit none
      
      character(len=*)        , intent(in)  :: param
      character(len=*)        , intent(in)  :: def
      logical                 , intent(out) :: isdef
      integer                               :: i
      
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Cannot get value. Reader not initialized! Assigning default value."
         getValue_st = def
         return
      end if
      
      getValue_st = def
      isdef       = .true.
      
      do i = 1, entries_found
         if (trim(ini_param(i)) .eq. trim(param)) then
            getValue_st = trim(ini_value(i))
            isdef       = .false.
            return
         end if
      end do
   end function getValue_st
         
   !****************************
   ! 
   !****************************  
   logical function getValue_lo(param, def, isdef)
      implicit none
      
      character(len=*), intent(in)  :: param
      logical         , intent(in)  :: def
      logical         , intent(out) :: isdef
      integer                       :: i
      
      if (.not. is_initialized) then
         write (*,*) "IniFileReader: Cannot get value. Reader not initialized! Assigning default value."
         getValue_lo = def
         return
      end if
      
      getValue_lo = def
      isdef       = .true.
      
      do i = 1, entries_found
         if (trim(ini_param(i)) .eq. trim(param)) then
            if ((trim(ini_value(i)) .eq. "true") .or. (trim(ini_value(i)) .eq. ".true.") .or. &
                (trim(ini_value(i)) .eq. "1") .or. (trim(ini_value(i)) .eq. "-1")) then
               getValue_lo = .true.
               isdef       = .false.
               return
            elseif ((trim(ini_value(i)) .eq. "false") .or. (trim(ini_value(i)) .eq. ".false.") .or. &
                    (trim(ini_value(i)) .eq. "0")) then
               getValue_lo = .false.
               isdef       = .false.
               return
            end if
         end if
      end do
   end function getValue_lo
end module ini_file_reader







