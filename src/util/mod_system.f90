!>
!! \brief Contains system-dependent procedures
!!
module mod_system
   implicit none

   character (len = *), parameter :: text_editor     = "emacs -nw "
   character (len = *), parameter :: graph_viewer    = "evince "
   character (len = *), parameter :: graph_generator = "gnuplot "

contains

   ! Gets the current process memory in MB
   ! This function was diactivated because it is not portable.
   ! In order to reactivate it, uncomment "write(PID,*) getPID()"
   real(8) function get_RAM(sim_id)
      implicit none
      character (len = *), intent(in) :: sim_id ! Simulation identification

      ! Inner variables

      character (len=20)  :: PID = ''
      character (len=100) :: fout
      real(8) :: RAM


      !write(PID,*) getPID()

      fout = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // "_memory.dat"

      call system("ps -o drs -p " // trim(adjustl(PID)) // " | tail -n 1 > " // trim(adjustl(fout)) )

      open(10,file=fout)

      read(10,*) RAM

      close(10)

      get_RAM = RAM / 1024.d0

   end function


   !> \brief Finds all the occurences of str1 in file1, replaces them by str2 and save the results in file2
   subroutine file_replace( file1, file2, str1, str2)
      implicit none
      character (len=*), intent(in) :: file1 !< Input file
      character (len=*), intent(in) :: file2 !< Output file
      character (len=*), intent(in) :: str1  !< Find
      character (len=*), intent(in) :: str2  !< Replace

      ! Auxiliary variables
      integer :: io
      character (len=1000) :: str

      if ( file1 /= file2 ) then

         open(10, file = file1)
         open(11, file = file2)

         do

            read(10,"(A)",iostat=io) str

            if ( io /= 0 ) exit

            call replace(str,str1,str2)

            write(11,*) trim(adjustl(str))

         end do

         close(10)
         close(11)

      else

         open(10, file = file1)
         open(11, file = trim(adjustl(file2))//".tmp")

         do

            read(10,"(A)",iostat=io) str

            if ( io /= 0 ) exit

            call replace(str,str1,str2)

            write(11,*) trim(adjustl(str))

         end do

         close(10)
         close(11)

         call rename( trim(adjustl(file2))//".tmp", file1 )

      end if

   end subroutine file_replace


   !> \brief Finds all the occurrences of str2 in str1 and replaces them by str3
   subroutine replace(str1, str2, str3)
      implicit none
      character (len=*), intent(inout) :: str1 !< Original string
      character (len=*), intent(in)    :: str2 !< Find
      character (len=*), intent(in)    :: str3 !< Replace

      ! Auxiliary variables

      integer :: idx !< First position where str2 was found
      integer :: n1  !< Length of str1
      integer :: n2  !< Length of str2

      n1 = len(str1)

      n2 = len(str2)

      do

         idx = index(str1,str2)

         if ( idx /= 0 ) then

            str1 = str1(1:idx-1) // str3 // str1(idx+n2:n1)

         else

            exit

         end if

      end do

   end subroutine replace

end module
