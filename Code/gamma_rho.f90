program calcular_gamma_rho

  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: max_points = 1000
  character(len=*), parameter :: default_output = "gamma_rho_points.dat"

  integer :: n
  integer :: i
  integer :: unit_in
  integer :: unit_out
  integer :: ios
  real(dp) :: phi_c
  real(dp) :: phi(max_points)
  real(dp) :: rho_c_max(max_points)
  real(dp) :: x(max_points)
  real(dp) :: y(max_points)
  real(dp) :: slope
  real(dp) :: gamma
  real(dp) :: intercept
  real(dp) :: r2
  character(len=256) :: input_file
  character(len=256) :: output_file
  character(len=512) :: line

  call read_arguments(phi_c, input_file, output_file)

  open(newunit=unit_in, file=trim(input_file), status="old", action="read", iostat=ios)
  if (ios /= 0) then
     write(*,*) "No pude abrir el archivo de entrada: ", trim(input_file)
     stop 1
  end if

  n = 0
  do
     read(unit_in,'(A)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0) cycle
     if (line(1:1) == "#") cycle

     n = n + 1
     if (n > max_points) error stop "Demasiados puntos en el archivo de entrada."

     read(line,*,iostat=ios) phi(n), rho_c_max(n)
     if (ios /= 0) then
        write(*,*) "Linea invalida en ", trim(input_file), ":"
        write(*,*) trim(line)
        stop 1
     end if

     if (phi(n) >= phi_c) then
        write(*,*) "Omito phi >= phi_c porque debe ser subcritico: ", phi(n)
        n = n - 1
        cycle
     end if
     if (rho_c_max(n) <= 0.0_dp) then
        write(*,*) "Omito rho_c_max no positiva para phi = ", phi(n)
        n = n - 1
        cycle
     end if

     x(n) = log(phi_c - phi(n))
     y(n) = log(rho_c_max(n))
  end do
  close(unit_in)

  if (n < 2) error stop "Se necesitan al menos dos puntos subcriticos validos."

  call linear_fit(n, x, y, slope, intercept, r2)
  gamma = -0.5_dp*slope

  open(newunit=unit_out, file=trim(output_file), status="replace", action="write")
  write(unit_out,'(A)') "# Ajuste: log(rho_c_max) = C - 2*gamma*log(phi_c - phi)"
  write(unit_out,'(A,ES24.16)') "# phi_c = ", phi_c
  write(unit_out,'(A,ES24.16)') "# pendiente = ", slope
  write(unit_out,'(A,ES24.16)') "# gamma = -pendiente/2 = ", gamma
  write(unit_out,'(A,ES24.16)') "# C = ", intercept
  write(unit_out,'(A,ES24.16)') "# R2 = ", r2
  write(unit_out,'(A)') "# phi rho_c_max log_phi_c_minus_phi log_rho_c_max"
  do i = 1, n
     write(unit_out,'(4ES24.15)') phi(i), rho_c_max(i), x(i), y(i)
  end do
  close(unit_out)

  write(*,'(A,I0)') "Puntos usados: ", n
  write(*,'(A,ES16.8)') "phi_c = ", phi_c
  write(*,'(A,ES16.8)') "pendiente = ", slope
  write(*,'(A,ES16.8)') "gamma = ", gamma
  write(*,'(A,ES16.8)') "C = ", intercept
  write(*,'(A,ES16.8)') "R2 = ", r2
  write(*,*) "Tabla guardada en: ", trim(output_file)

contains

  subroutine read_arguments(phi_c, input_file, output_file)

    real(dp), intent(out) :: phi_c
    character(len=*), intent(out) :: input_file
    character(len=*), intent(out) :: output_file
    character(len=256) :: arg
    integer :: nargs

    nargs = command_argument_count()
    if (nargs < 2 .or. nargs > 3) then
       write(*,*) "Uso:"
       write(*,*) "  calcular_gamma_rho.exe phi_c rho_runs.dat [gamma_rho_points.dat]"
       write(*,*)
       write(*,*) "Formato de rho_runs.dat:"
       write(*,*) "  # phi rho_c_max"
       write(*,*) "  0.352900 1.2345e6"
       write(*,*)
       write(*,*) "Usa solamente corridas subcriticas: phi < phi_c."
       stop
    end if

    call get_command_argument(1, arg)
    read(arg,*) phi_c
    call get_command_argument(2, input_file)

    if (nargs == 3) then
       call get_command_argument(3, output_file)
    else
       output_file = default_output
    end if

  end subroutine read_arguments

  subroutine linear_fit(n, x, y, slope, intercept, r2)

    integer, intent(in) :: n
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: slope
    real(dp), intent(out) :: intercept
    real(dp), intent(out) :: r2
    real(dp) :: sx
    real(dp) :: sy
    real(dp) :: sxx
    real(dp) :: sxy
    real(dp) :: ymean
    real(dp) :: ss_tot
    real(dp) :: ss_res
    real(dp) :: denom
    real(dp) :: yfit
    integer :: i

    sx = sum(x(1:n))
    sy = sum(y(1:n))
    sxx = sum(x(1:n)*x(1:n))
    sxy = sum(x(1:n)*y(1:n))
    denom = real(n, dp)*sxx - sx*sx
    if (abs(denom) <= tiny(1.0_dp)) error stop "No se puede ajustar: los puntos tienen la misma x."

    slope = (real(n, dp)*sxy - sx*sy)/denom
    intercept = (sy - slope*sx)/real(n, dp)

    ymean = sy/real(n, dp)
    ss_tot = 0.0_dp
    ss_res = 0.0_dp
    do i = 1, n
       yfit = intercept + slope*x(i)
       ss_tot = ss_tot + (y(i) - ymean)**2
       ss_res = ss_res + (y(i) - yfit)**2
    end do

    if (ss_tot > 0.0_dp) then
       r2 = 1.0_dp - ss_res/ss_tot
    else
       r2 = 1.0_dp
    end if

  end subroutine linear_fit

end program calcular_gamma_rho
