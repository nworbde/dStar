module NScool_epochs
    use NScool_def
    integer, parameter :: initial_epoch_length = 64
    integer, parameter :: token_null   = 0
    integer, parameter :: token_EOF    = 1
    integer, parameter :: token_string = 2
    integer, parameter :: token_value  = 3

contains
    
    subroutine do_load_epochs(s, ierr)
        use iso_fortran_env, only : IOSTAT_END
        use exceptions_lib
        use constants_def, only : julian_day
        use storage, only: allocate_NScool_epoch_arrays
        
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine = 'do_load_epochs'
        character(len=256) :: buffer, string, bad_line_msg
        integer :: file_id, t, n, i, m, i0, i1
        integer :: lines_read,iostat,this_line,epoch_length,n_cycles,bad_lines
        integer :: col_time, col_Mdot
        real(dp) :: delta_t
        real(dp), dimension(:,:), allocatable :: in_array
        type(failure) :: file_error=failure(scope='do_load_epochs')
        type(alert) :: file_not_empty=alert(&
        &   scope='do_load_epochs',message='unexpected EOF')
        type(alert) :: bad_token=alert(scope='do_load_epochs')
        type(failure) :: allocation_error=failure(scope='do_load_epochs', &
        &   message='allocating storage')

        if (.not. s% load_epochs) return

        ! set defaults
        ierr = 0
        n = 0
        i = 0
        n_cycles = 1
        col_time = 1
        col_Mdot = 2

        open(newunit=file_id,file=trim(s% epoch_datafile), &
            & status='old',action='read',iostat=ierr)
        if (file_error% raised(ierr, &
        &   'unable to open ' // trim(s% epoch_datafile))) return
        
        ! there are 5 data flags
        !   number_cycles: allows for a repeating pattern
        !   Mdot_scale: useful if the table is in Eddington values and need to 
        !       be converted to cgs scale
        !   time_scale: for converting time to solar days
        !   columns: indicates order of columns, default is time Mdot
        !   epochs: loader assumes rest of file are columns of time [MJD],
        !       accretion rate
        scan_file: do
            t = token(file_id,n,i,buffer,string)
            if (t == token_EOF) then
                ! note this does not raise an error flag, but the number of
                ! epochs will be 0
                call file_not_empty% report
                return
            end if
            if (t == token_value) then
                select case(trim(string))
                case('number_cycles')
                    t = token(file_id,n,i,buffer,string,advance=.FALSE.)
                    if (t /= token_value) then
                        call bad_token% report('reading number_cycles')
                        cycle
                    end if
                    read(string,*) n_cycles
                case('Mdot_scale')
                    t = token(file_id,n,i,buffer,string,advance=.FALSE.)
                    if (t /= token_value) then
                        call bad_token% report('reading Mdot_scale')
                        cycle
                    end if
                    read(string,*) s% epoch_Mdot_scale
                case('time_scale')
                    t = token(file_id,n,i,buffer,string,advance=.FALSE.)
                    if (t /= token_value) then
                        call bad_token% report('reading time_scale')
                        cycle
                    end if
                    read(string,*) s% epoch_time_scale
                case('columns')
                    t = token(file_id,n,i,buffer,string,advance=.FALSE.)
                    if (t /= token_string) then
                        call bad_token% report('reading columns')
                        cycle
                    end if
                    if (string(1:4) == 'Mdot') then
                        col_Mdot = 1; col_time = 2
                    else if (string(1:4) == 'time') then
                        col_time = 1; col_Mdot = 2
                    else
                        call bad_token% report('reading columns')
                        cycle
                    end if
                case('epochs')
                    exit    ! terminate cycle, read in data table
                case default
                    call bad_token% report('unexpected token ' // trim(string))
                end select
            end if
        end do scan_file
        
        epoch_length = initial_epoch_length
        allocate(in_array(2,epoch_length), stat=ierr)
        if (allocation_error% raised(ierr)) return
        
        lines_read = 0
        bad_lines = 0
        read_table: do
            this_line = lines_read+1
            read(file_id,*,iostat=ierr) in_array(:,this_line)
            if (ierr == IOSTAT_END) then
                ierr = 0
                exit
            else if (ierr /= 0) then
                ! bad line, issue warning and skip
                write(bad_line_msg,'(a,tr1,i0,tr1,a)') 'unable to read line',this_line+bad_lines,'of data; skipping'
                call bad_token% report(bad_line_msg)
                bad_lines = bad_lines+1
                cycle
            end if
            lines_read = this_line
            if (this_line == size(in_array,dim=2)) then
                epoch_length = 2*size(in_array,dim=2)
                call realloc(in_array,2,epoch_length,ierr)
                if (allocation_error% raised(ierr)) exit
            end if
        end do read_table
        close(file_id)
        
        ! scale the accretion rates and times
        in_array(col_time,:) = in_array(col_time,:) * &
        &   s% epoch_time_scale * julian_day
        in_array(col_Mdot,:) = in_array(col_Mdot,:) * s% epoch_Mdot_scale
        
        if (ierr == 0) then
            ! Store the table. For n lines read and m cycles, we repeat
            ! lines 1:n-1 m times. The times are stored in epoch_boundares 
            ! 0:(n-1)-1, n-1:2(n-1)-1, ..., (m-1)*(n-1):m(n-1)-1;
            ! if the table has times t_1,...,t_n, then with Delta = t_n-t_1 the 
            ! stored times are
            ! t_1,t_2,...,t_(n-1),t_1+Delta,...,t_(n-1)+Delta,t_1+2*Delta,...
            ! t_(n-1)+(m-1)*Delta. The final time is stored in 
            ! epoch_boundaries(m*(n-1)) = t_n+(m-1)*Delta.
            ! The accretion rates are stored in
            ! epoch_Mdots 1:n-1, (n-1)+1:2(n-1), ..., (m-1)*(n-1)+1:m(n-1).
            ! 
            n = lines_read
            s% number_epochs = n_cycles*(n-1)
            call allocate_NScool_epoch_arrays(s, ierr)
            delta_t = in_array(col_time,n)-in_array(col_time,1)
            do m = 0, n_cycles-1
                i0 = m*(n-1)
                i1 = i0+n-2
                s% epoch_boundaries(i0:i1) = in_array(col_time,1:n-1)  &
                &   + m*delta_t
                s% epoch_Mdots(i0+1:i1+1) =  in_array(col_Mdot,1:n-1)
            end do
            ! final time (end of integration)
            s% epoch_boundaries(s% number_epochs) = in_array(col_time,n) &
            &   + (n_cycles-1)*delta_t
        end if
        
        deallocate(in_array)
     end subroutine do_load_epochs
     
     subroutine realloc(a,n1,n2,ierr)
         real(dp), dimension(:,:), allocatable, intent(inout) :: a
         integer, intent(in) :: n1, n2
         integer, intent(out) :: ierr
         real(dp), dimension(size(a,dim=1),size(a,dim=2)) :: tmp
         integer :: m1, m2
         
         m1 = size(a,dim=1)
         m2 = size(a,dim=2)
         ierr = 0
         tmp = a
         deallocate(a)
         allocate(a(n1,n2),stat = ierr)
         if (ierr /= 0) return
         a(1:m1,1:m2) = tmp(1:m1,1:m2)
     end subroutine realloc

     function token(fid,n,i,buffer,string,advance)
         use iso_fortran_env, only : IOSTAT_END
         integer, intent(in) :: fid
         integer, intent(inout) :: n,i
         character(len=*), intent(inout) :: buffer,string
         logical, intent(in), optional :: advance
         integer :: token
         logical :: okay_to_advance
         integer :: info, j, j1, j2, l, str_len
         character(len=1), parameter :: tab = char(9)
         character(len=1) :: delim

         okay_to_advance = .TRUE.
         if (present(advance)) okay_to_advance = advance
         
         token = token_null
         line_loop: do
             do while (i >= n)
                 if (.not.okay_to_advance) return
                 read(fid,fmt='(a)',iostat=info) buffer
                 if (info == IOSTAT_END) then
                     token = token_EOF
                     return
                 end if
                 n = len_trim(buffer)
                 i = 0
             end do
             
             token_loop: do while (i < n)
                 i = i + 1
                 select case(buffer(i:i))
                 case ('!')
                     i = n
                     cycle line_loop
                 case (' ',tab)
                     cycle token_loop
                 case ('"','''')
                     j = 1; str_len = len(string)
                     delim = buffer(i:i)
                     string_loop: do
                         i = i + 1
                         if (i > n) exit
                         if (buffer(i:i) == delim) exit
                         if (j > str_len) exit
                         string(j:j) = buffer(i:i)
                         j = j+1
                     end do string_loop
                     do while (j <= str_len)
                        string(j:j) = ' '
                        j = j+1
                     end do
                     token = token_string
                     return
                 case default
                     j1 = i; j2 = i; str_len = len(string)
                     value_loop: do
                         if (i+1 > n) exit
                         if (buffer(i+1:i+1) == ' ') exit
                         if (buffer(i+1:i+1) == tab) exit
                         i = i + 1
                         j2 = i
                     end do value_loop
                     l = min(j2-j1+1,str_len)
                     string(1:l) = buffer(j1:j2)
                     do j = l+1,str_len
                         string(j:j) = ' '
                     end do
                     token = token_value
                     return
                 end select
             end do token_loop
         end do line_loop
     end function token
end module NScool_epochs
