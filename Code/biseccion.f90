program biseccion_phi0

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  character(len=*), parameter :: collapse_lapse_text = &
       "Resultado: Colapso de la funcion lapso."
  character(len=*), parameter :: disperse_lapse_text = &
       "Resultado: evolucion terminada sin colapso de la funcion lapso."
  character(len=*), parameter :: collapse_geometry_text = &
       "Resultado: Colapso geometrico por max(2m/r)."
  character(len=*), parameter :: disperse_geometry_text = &
       "Resultado: evolucion terminada sin colapso geometrico."

  integer :: nr_arg
  integer :: max_iter
  integer :: it
  integer :: narg
  real(dp) :: rmax
  real(dp) :: t_final
  real(dp) :: weak
  real(dp) :: strong
  real(dp) :: mid
  real(dp) :: tol
  real(dp) :: delta
  real(dp) :: estimate
  character(len=32) :: status
  character(len=256) :: out_dir
  character(len=256) :: temp_log_file
  character(len=256) :: csv_file

  narg = command_argument_count()
  if (narg < 5 .or. narg > 7) then
     call print_usage()
     stop
  end if

  call read_integer_argument(1, nr_arg)
  call read_real_argument(2, rmax)
  call read_real_argument(3, t_final)
  call read_real_argument(4, weak)
  call read_real_argument(5, strong)

  tol = 1.0d-14
  if (narg >= 6) call read_real_argument(6, tol)

  max_iter = 80
  if (narg >= 7) call read_integer_argument(7, max_iter)

  if (weak >= strong) error stop "El valor weak debe ser menor que strong."
  if (tol <= 0.0_dp) error stop "La tolerancia debe ser positiva."
  if (max_iter < 1) error stop "max_iter debe ser positivo."

  out_dir = "biseccion_phi0_fortran"
  csv_file = trim(out_dir)//"\biseccion.csv"
  temp_log_file = trim(out_dir)//"\run_actual.tmp"

  call make_directory(out_dir)

  write(*,*) "Validando extremos iniciales..."

  call run_and_classify(nr_arg, rmax, t_final, weak, temp_log_file, status)
  if (trim(status) /= "disperse") then
     write(*,*) "El extremo weak no disperso. Resultado: ", trim(status)
     error stop
  end if

  call run_and_classify(nr_arg, rmax, t_final, strong, temp_log_file, status)
  if (trim(status) /= "collapse") then
     write(*,*) "El extremo strong no colapso. Resultado: ", trim(status)
     error stop
  end if

  open(unit=20, file=trim(csv_file), status="replace", action="write")
  write(20,'(A)') "iter,weak,strong,mid,resultado,delta_absoluto,delta_relativo"

  do it = 1, max_iter
     mid = 0.5_dp*(weak + strong)

     if (mid <= weak .or. mid >= strong) then
        delta = (strong - weak)/abs(weak)
        write(*,*) "El intervalo ya no se puede partir con precision double."
        exit
     end if

     call run_and_classify(nr_arg, rmax, t_final, mid, temp_log_file, status)

     select case (trim(status))
     case ("disperse")
        weak = mid
     case ("collapse")
        strong = mid
     case default
        delta = (strong - weak)/abs(weak)
        write(20,'(I0,",",ES26.18,",",ES26.18,",",ES26.18,",",A,",",ES26.18,",",ES26.18)') &
             it, weak, strong, mid, trim(status), strong - weak, delta
        close(20)
        write(*,*) "No pude clasificar la iteracion ", it
        error stop
     end select

     delta = (strong - weak)/abs(weak)
     estimate = 0.5_dp*(weak + strong)

     write(20,'(I0,",",ES26.18,",",ES26.18,",",ES26.18,",",A,",",ES26.18,",",ES26.18)') &
          it, weak, strong, mid, trim(status), strong - weak, delta
     flush(20)

     write(*,'(A,I0,A,ES26.18,A,A,A,ES26.18,A,ES10.3)') &
          "iter=", it, " mid=", mid, " resultado=", trim(status), &
          " phi_crit~", estimate, " delta=", delta

     if (delta < tol) exit
  end do

  close(20)

  write(*,*)
  write(*,*) "Tabla guardada en: ", trim(csv_file)
  write(*,'(A,ES26.18,A,ES26.18,A)') "Intervalo final: [", weak, ", ", strong, "]"
  write(*,'(A,ES26.18)') "phi0 critico estimado: ", 0.5_dp*(weak + strong)
  write(*,'(A,ES10.3)') "delta relativo: ", (strong - weak)/abs(weak)

contains

  subroutine print_usage()
    write(*,*) "Uso:"
    write(*,*) "  biseccion_phi0.exe Nr_intervalos rmax t_final weak strong [tol] [max_iter]"
    write(*,*) "Ejemplo:"
    write(*,*) "  biseccion_phi0.exe 100 10.0 6.0 0.40 0.50 1e-14 80"
    write(*,*)
    write(*,*) "weak  = valor que NO colapsa"
    write(*,*) "strong = valor que SI colapsa"
  end subroutine print_usage

  subroutine read_integer_argument(iarg, value)
    integer, intent(in) :: iarg
    integer, intent(out) :: value
    character(len=128) :: arg

    call get_command_argument(iarg, arg)
    read(arg,*) value
  end subroutine read_integer_argument

  subroutine read_real_argument(iarg, value)
    integer, intent(in) :: iarg
    real(dp), intent(out) :: value
    character(len=128) :: arg

    call get_command_argument(iarg, arg)
    read(arg,*) value
  end subroutine read_real_argument

  subroutine make_directory(path)
    character(len=*), intent(in) :: path
    character(len=512) :: command
    integer :: exitstat
    integer :: cmdstat

    command = 'cmd /c if not exist "'//trim(path)//'" mkdir "'//trim(path)//'"'
    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)

    if (cmdstat /= 0 .or. exitstat /= 0) then
       write(*,*) "No se pudo crear el directorio: ", trim(path)
       error stop
    end if
  end subroutine make_directory

  subroutine run_and_classify(nr_arg, rmax, t_final, phi0, log_file, status)
    integer, intent(in) :: nr_arg
    real(dp), intent(in) :: rmax
    real(dp), intent(in) :: t_final
    real(dp), intent(in) :: phi0
    character(len=*), intent(in) :: log_file
    character(len=*), intent(out) :: status

    character(len=1024) :: command
    integer :: exitstat
    integer :: cmdstat

    write(command,'(A,I0,A,ES24.16,A,ES24.16,A,ES24.16,A,A,A)') &
         'cmd /c .\CEesferico1D.exe ', nr_arg, ' ', rmax, ' ', t_final, ' ', phi0, &
         ' > "', trim(log_file), '" 2>&1'

    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)

    if (cmdstat /= 0 .or. exitstat /= 0) then
       call delete_file(log_file)
       status = "error"
       return
    end if

    call classify_log(log_file, status)
    call delete_file(log_file)
  end subroutine run_and_classify

  subroutine classify_log(log_file, status)
    character(len=*), intent(in) :: log_file
    character(len=*), intent(out) :: status

    character(len=512) :: line
    integer :: ios
    logical :: found_collapse
    logical :: found_disperse

    found_collapse = .false.
    found_disperse = .false.

    open(unit=30, file=trim(log_file), status="old", action="read", iostat=ios)
    if (ios /= 0) then
       status = "error"
       return
    end if

    do
       read(30,'(A)', iostat=ios) line
       if (ios /= 0) exit
       if (index(line, collapse_lapse_text) > 0 .or. &
            index(line, collapse_geometry_text) > 0) found_collapse = .true.
       if (index(line, disperse_lapse_text) > 0 .or. &
            index(line, disperse_geometry_text) > 0) found_disperse = .true.
    end do

    close(30)

    if (found_collapse) then
       status = "collapse"
    else if (found_disperse) then
       status = "disperse"
    else
       status = "unknown"
    end if
  end subroutine classify_log

  subroutine delete_file(path)
    character(len=*), intent(in) :: path
    character(len=512) :: command
    integer :: exitstat
    integer :: cmdstat

    command = 'cmd /c if exist "'//trim(path)//'" del /q "'//trim(path)//'"'
    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)
  end subroutine delete_file

end program biseccion_phi0
