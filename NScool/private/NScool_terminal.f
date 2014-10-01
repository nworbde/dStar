module NScool_terminal
   use NScool_def
   
   character(len=*), parameter :: scival = 'es15.6', intval = 'i15', fltval = 'f15.6'
   integer, parameter :: num_terminal_cols = 11, terminal_col_width=15

   character(len=*), parameter :: terminal_count_fmt = '(10'//intval//')'
   character(len=*), parameter :: terminal_title_fmt = '(10a15)'
   character(len=*), parameter :: terminal_val_fmt = '('//intval//',2'//scival//',8'//fltval//')'
   character(len=terminal_col_width),dimension(num_terminal_cols) :: terminal_cols  &
         & = [character(len=terminal_col_width) ::  'model', 'time', 'dt', 'lg(Lsurf)', 'lg(Teff)', &
         & 'lg(Lnu)','lg(Lnuc)','max(lgT)','lg(P[Tmax])','min(lgT)','lg(P[Tmin])' ]
   
   ! logical :: terminal_first_call = .TRUE.
   
   contains
   subroutine do_write_terminal(id,ierr,print_header)
       use constants_def
      use iso_fortran_env, only : output_unit
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      logical, intent(in), optional :: print_header
      type(NScool_info), pointer :: s
      real :: Lnu, Lnuc
      character(len=16) ::terminal_date, terminal_time, terminal_zone

      call get_NScool_info_ptr(id,s,ierr)
      if (ierr /= 0) return
            
      Lnu = dot_product(s% enu,s% dm)
      Lnuc = dot_product(s% enuc, s% dm)

      ! got rid of the datestring write, for purposes of testing output
      ! first time called: write datestring
      ! if (terminal_first_call) then
      !     call date_and_time(terminal_date,terminal_time,terminal_zone)
      !     write (output_unit,'(a,"/",2a,/)') trim(terminal_date),trim(terminal_time),trim(terminal_zone)
      !     terminal_first_call = .FALSE.
      !  end if
       
      ! write a header line upon request
      if (present(print_header)) then
         if (print_header) then
            write (output_unit,*)
            write(output_unit,terminal_title_fmt) adjustr(terminal_cols)
            write (output_unit,*)
         end if
      end if
      
      ! write the row
      write(output_unit,terminal_val_fmt) s% model, s% tsec, s% dt, log10(s% Lsurf), log10(s% Teff),  &
         &  log10(Lnu), log10(Lnuc), maxval(s% lnT)/ln10, log10(s% P(maxloc(s% lnT))),  &
         &  minval(s% lnT)/ln10, log10(s% P(minloc(s% T)))
      
   end subroutine do_write_terminal

end module NScool_terminal

