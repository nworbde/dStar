module NScool_private_def
   use NScool_def

   character(len=256) :: dStar_data_dir
   
   contains
   
   subroutine NScool_private_def_init(my_dStar_dir)
      character(len=*), intent(in) :: my_dStar_dir 
      integer :: i
      do i=1,max_NScool_handles
         NScool_handles(i)% id = i
         NScool_handles(i)% in_use = .false.
      end do
      
      dStar_data_dir = trim(my_dStar_dir)//'/data'
   end subroutine NScool_private_def_init

   subroutine free_NScool(s)
      type(NScool_info), pointer :: s
      NScool_handles(s% id)% in_use = .FALSE.
   end subroutine free_NScool

   function NScool_private_alloc(ierr)
      integer, intent(out) :: ierr
      integer :: NScool_private_alloc
      integer :: i
      
      ierr = 0
      NScool_private_alloc = -1
      do i=1, max_NScool_handles
         if (.not. NScool_handles(i)% in_use) then
            NScool_handles(i)% in_use = .TRUE.
            NScool_private_alloc = i
            exit
         end if
      end do
      if (NScool_private_alloc == -1) then
         ierr = -1
         return
      end if
      if (NScool_handles(NScool_private_alloc)% id /= NScool_private_alloc) then
         ierr = -1
         return
      end if
   end function NScool_private_alloc

end module NScool_private_def
