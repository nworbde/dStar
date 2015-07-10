module NScool_history
   use NScool_def
   
   character(len=*), parameter :: scival = 'es15.6', intval = 'i15', fltval = 'f15.6'
   integer, parameter :: num_header_cols = 4, header_col_width=15
   integer, parameter :: num_history_cols = 7, history_col_width=15

   character(len=*), parameter :: header_count_fmt = '(4'//intval//')'
   character(len=*), parameter :: header_title_fmt = '(4a15)'
   character(len=*), parameter :: header_val_fmt = '(4'//scival//')'
   character(len=*), parameter :: history_count_fmt = '(7'//intval//')'
   character(len=*), parameter :: history_title_fmt = '(7a15)'
   character(len=*), parameter :: history_val_fmt = '('//intval//scival//',5'//fltval//')'
   character(len=header_col_width), dimension(num_header_cols) :: header_cols = [character(len=header_col_width) ::  &
      & 'gravity', 'core_mass', 'core_radius', 'core_temp.' ]
   character(len=history_col_width),dimension(num_history_cols) :: history_cols = [character(len=history_col_width) ::  &
      & 'model','time','lg(Mdot)','lg(Teff)','lg(Lsurf)','lg(Lnu)','lg(Lnuc)' ]
   
   logical, save :: history_first_call = .TRUE.
   
   contains
   subroutine do_write_history(id,ierr)
      use utils_lib, only: alloc_iounit, free_iounit, append_line
      use constants_def, only: julian_day
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(NScool_info), pointer :: s
      integer :: iounit, iz, ih
      real(dp) :: Lnu, Lnuc
      character(len=16) ::history_date, history_time, history_zone
      
      call get_NScool_info_ptr(id,s,ierr)
      if (ierr /= 0) return
      
      if (history_first_call) then
      
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
      
         open(unit=iounit,file=trim(s% history_filename),action='write',iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open history file ',trim(s% history_filename)
            call free_iounit(iounit)
            return
         end if
               
         ! first line: write datestring
         call date_and_time(history_date,history_time,history_zone)
         write (iounit,'(a,"/",2a,/)') trim(history_date),trim(history_time),trim(history_zone)
         write (iounit,header_count_fmt) (ih, ih=1,num_header_cols)
         write (iounit,header_title_fmt) adjustr(header_cols)
         write (iounit,header_val_fmt) s% grav, s% Mcore, s% Rcore, s% Tcore
         write (iounit,*)
         write (iounit,*)
         write(iounit,history_count_fmt) (ih, ih=1,num_history_cols)
         write(iounit,history_title_fmt) adjustr(history_cols)
         close(iounit)
         call free_iounit(iounit)
         history_first_call = .FALSE.
         
      end if
         
      ! subsequent calls
      Lnu = dot_product(s% enu, s% dm)
      Lnuc = dot_product(s% enuc, s% dm)
        
      iounit = alloc_iounit(ierr)
      if (ierr /= 0) return
   
      open(unit=iounit,file=trim(s% history_filename),action='write',position='append',iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'failed to open history file ',trim(s% history_filename)
         call free_iounit(iounit)
         return
      end if
      write(iounit,history_val_fmt) s% model, s% tsec/julian_day, safelog10(s% Mdot), log10(s% Teff), log10(s% Lsurf),  &
      & log10(Lnu), safelog10(Lnuc)
      close(iounit)
      call free_iounit(iounit)
  contains
      function safelog10(x)
          real(dp), intent(in) :: x
          real(dp) :: safelog10
          if (x <= 0.0_dp) then
              safelog10 = 0.0_dp
          else 
              safelog10 = log10(x)
          end if
      end function safelog10
   end subroutine do_write_history

end module NScool_history

