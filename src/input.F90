! ##################################
! INPUT ROUTINES
! ##################################
module input
   use global
   use output_data
   use ini_file_reader
   
   implicit none
   
   interface getInputValue        
      module procedure getInputValue_sp     ! single precision
      module procedure getInputValue_dp     ! double precision
      module procedure getInputValue_i      ! integer
      module procedure getInputValue_st     ! string
      module procedure getInputValue_lo     ! logical
   end interface
   
   ! these are the variables that can be passed on the command line
   integer :: cl_ompThreads
   
   contains
   
   !------------------------------
   ! Parse the command line arguments
   !------------------------------   
   subroutine parseCommandLine()
      implicit none
      
      integer                           :: i, nargs
      character(len=200)                :: commandLineArgument
      character(len=200), dimension(20) :: args
      
      ! first we have to initialize some variables to a default so we can check
      ! when we parse the input file if the variable has been given on the command line
      cl_ompThreads = 0
      
      
      ! first argument has to be the input file
      call get_command_argument(1, inputFile)
      
      ! check for additional arguments
      do i = 2, command_argument_count()
         call get_command_argument(i, commandLineArgument)
         
         ! split the command line argument at "=" in the array args
         call parse(commandLineArgument, "=", args, nargs)
         
         select case (trim(args(1)))
            case ("ncpus")
               read (args(2), *) cl_ompThreads
            case default
               write (*,'(2A)') " WARNING: Unrecognized command-line argument: ", trim(args(1))
         end select
      end do
   
   end subroutine parseCommandLine
   
   
   !------------------------------
   ! Read the input file
   !------------------------------   
   subroutine readInputFile(inputFile)
      implicit none
      
      character(len=*), intent(in) :: inputFile
      integer                  :: ios
      logical                  :: fileExists
      
      ! check & open the input file
      inquire(file=trim(inputFile), exist=fileExists, iostat=ios)

      if (.not. fileExists) then
         write (*,*) "Input file ", trim(inputFile), " not found! STOP"
         !call MPI_FINALIZE (mpiErr) 
         stop
      end if   
      
      open (unit=fid_Input, action="read", file=trim(inputFile), iostat=ios)

      if (ios .ne. 0) then
         write (*,*) "Could not open the input file! STOP"
         !call MPI_FINALIZE (mpiErr) 
         stop
      end if
      
      ! analyze the input file
      call initIniFileReader(fid_Input, 300)
      call analyzeIniFile(ios)
      
      if (ios .ne. 0) then
         write (*,*) "There was a problem with the input file! STOP"
         call finishIniFileReader()
         close (fid_Input)
         !call MPI_FINALIZE (mpiErr) 
         stop
      end if
      
      close (fid_Input)
      
   end subroutine readInputFile
      
   !------------------------------
   ! Get a parameter from the input file and echo it
   !------------------------------   
   real function getInputValue_sp(param, def)
      implicit none
      
      character(len=*), intent(in)  :: param
      real            , intent(in)  :: def
      logical                       :: isdef
      real                          :: val
      
      val = getValue(param, def, isdef)
      
      write (*, '(A7,A35,A3,E16.7)', advance='no') "Input: ", trim(param), " = ", val
      if (isdef) then
         write (*,'(A)') , "        DEFAULT!"
      else
         write (*,*)
      end if     
      
      getInputValue_sp = val
      
   end function getInputValue_sp   
   
   !------------------------------
   ! Get a parameter from the input file and echo it
   !------------------------------   
   double precision function getInputValue_dp(param, def)
      implicit none
      
      character(len=*), intent(in)  :: param
      double precision, intent(in)  :: def
      logical                       :: isdef
      double precision              :: val
      
      val = getValue(param, def, isdef)
      
      write (*, '(A7,A35,A3,E16.7)', advance='no') "Input: ", trim(param), " = ", val
      if (isdef) then
         write (*,'(A)') , "        DEFAULT!"
      else
         write (*,*)
      end if     
      
      getInputValue_dp = val
      
   end function getInputValue_dp
      
   !------------------------------
   ! Get a parameter from the input file and echo it
   !------------------------------   
   integer function getInputValue_i(param, def)
      implicit none
      
      character(len=*), intent(in)  :: param
      integer         , intent(in)  :: def
      logical                       :: isdef
      integer                       :: val
      
      val = getValue(param, def, isdef)
      
      write (*, '(A7,A35,A3,I16)', advance='no') "Input: ", trim(param), " = ", val
      if (isdef) then
         write (*,'(A)') , "        DEFAULT!"
      else
         write (*,*)
      end if     
      
      getInputValue_i = val
      
   end function getInputValue_i
      
   !------------------------------
   ! Get a parameter from the input file and echo it
   !------------------------------   
   character(len=value_len) function getInputValue_st(param, def)
      implicit none
      
      character(len=*), intent(in)  :: param
      character(len=*), intent(in)  :: def
      logical                       :: isdef
      character(len=value_len)      :: val
      
      val = getValue(param, def, isdef)
      
      write (*, '(A7,A35,A3,A16)', advance='no') "Input: ", trim(param), " = ", trim(val)
      if (isdef) then
         write (*,'(A)') , "        DEFAULT!"
      else
         write (*,*)
      end if     
      
      getInputValue_st = trim(val)
      
   end function getInputValue_st   
      
   !------------------------------
   ! Get a parameter from the input file and echo it
   !------------------------------   
   logical function getInputValue_lo(param, def)
      implicit none
      
      character(len=*), intent(in)  :: param
      logical         , intent(in)  :: def
      logical                       :: isdef
      logical                       :: val
      
      val = getValue(param, def, isdef)
      
      write (*, '(A7,A35,A3,L16)', advance='no') "Input: ", trim(param), " = ", val
      if (isdef) then
         write (*,'(A)') , "        DEFAULT!"
      else
         write (*,*)
      end if     
      
      getInputValue_lo = val
      
   end function getInputValue_lo     
end module input
