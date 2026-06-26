program CEesferico1D

  use evolution
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: OUTPUT_EVERY = 1
  real(dp), parameter :: CFL_FACTOR = 0.50_dp
  real(dp), parameter :: LAPSE_COLLAPSE_CRITERION = 1.0e-3_dp
  real(dp), parameter :: r0 = 0.0_dp
  real(dp), parameter :: sigma = 1.0_dp

  integer :: Nr
  integer :: Nt
  integer :: j
  integer :: n
  real(dp) :: rmax
  real(dp) :: t_final
  real(dp) :: phi0
  character(len=100) :: phi0_input
  integer :: phi0_decimals
  real(dp) :: dr
  real(dp) :: dt
  real(dp) :: t
  real(dp) :: rho_c_max
  real(dp) :: rho_c_max_time
  logical :: lapse_collapsed

  real(dp), allocatable :: r(:)
  real(dp), allocatable :: scalar_nm1(:), scalar_n(:), scalar_np1(:)
  real(dp), allocatable :: Phi_nm1(:), Phi_n(:), Phi_np1(:)
  real(dp), allocatable :: Pi_nm1(:), Pi_n(:), Pi_np1(:)
  real(dp), allocatable :: a_nm1(:), a_n(:), a_np1(:)
  real(dp), allocatable :: alpha_nm1(:), alpha_n(:), alpha_np1(:)

  call read_parameters(Nr, rmax, t_final, phi0, phi0_input, phi0_decimals)
  call delete_previous_outputs()

  dr = rmax/real(Nr - 1, dp)
  dt = CFL_FACTOR*dr
  Nt = int(t_final/dt + 0.5_dp)

  allocate(r(0:Nr-1))
  allocate(scalar_nm1(0:Nr-1), scalar_n(0:Nr-1), scalar_np1(0:Nr-1))
  allocate(Phi_nm1(0:Nr-1), Phi_n(0:Nr-1), Phi_np1(0:Nr-1))
  allocate(Pi_nm1(0:Nr-1), Pi_n(0:Nr-1), Pi_np1(0:Nr-1))
  allocate(a_nm1(0:Nr-1), a_n(0:Nr-1), a_np1(0:Nr-1))
  allocate(alpha_nm1(0:Nr-1), alpha_n(0:Nr-1), alpha_np1(0:Nr-1))

  do j = 0, Nr - 1
     r(j) = real(j, dp)*dr
  end do

  scalar_nm1 = 0.0_dp
  scalar_n = 0.0_dp
  scalar_np1 = 0.0_dp
  Phi_nm1 = 0.0_dp
  Phi_n = 0.0_dp
  Phi_np1 = 0.0_dp
  Pi_nm1 = 0.0_dp
  Pi_n = 0.0_dp
  Pi_np1 = 0.0_dp
  a_nm1 = 0.0_dp
  a_n = 0.0_dp
  a_np1 = 0.0_dp
  alpha_nm1 = 0.0_dp
  alpha_n = 0.0_dp
  alpha_np1 = 0.0_dp

  call initial_condition(Nr, r, dr, phi0, r0, sigma, &
       scalar_nm1, Phi_nm1, Pi_nm1, a_nm1, alpha_nm1)

  t = 0.0_dp
  rho_c_max = central_energy_density(Phi_nm1, Pi_nm1, a_nm1)
  rho_c_max_time = t
  lapse_collapsed = .false.

  call output_snapshot(0, Nr, r, scalar_nm1, Phi_nm1, Pi_nm1, a_nm1, alpha_nm1)

  Phi_n(0) = 0.0_dp
  a_n(0) = 1.0_dp
  alpha_n(0) = 1.0_dp

  call time_step(0, Nr, r, dr, dt, &
       scalar_nm1, Phi_nm1, Pi_nm1, a_nm1, alpha_nm1, &
       Phi_nm1, Pi_nm1, a_nm1, alpha_nm1, &
       Phi_n, Pi_n, scalar_n)

  call radiation_boundary_condition(0, Nr, r, dr, dt, &
       scalar_nm1, Pi_nm1, scalar_nm1, Phi_nm1, a_nm1, alpha_nm1, &
       scalar_n, Phi_n, Pi_n)

  Pi_n(0) = -Pi_nm1(0) + Pi_n(1) + Pi_nm1(1)
  call solve_metric(Nr, r, dr, Phi_n, Pi_n, a_n, alpha_n)
  t = t + 0.5_dp*dt
  call update_central_density_max(t, Phi_n, Pi_n, a_n, rho_c_max, rho_c_max_time)

  Phi_np1(0) = 0.0_dp
  a_np1(0) = 1.0_dp
  alpha_np1(0) = 1.0_dp

  call time_step(1, Nr, r, dr, dt, &
       scalar_n, Phi_n, Pi_n, a_n, alpha_n, &
       Phi_nm1, Pi_nm1, a_nm1, alpha_nm1, &
       Phi_np1, Pi_np1, scalar_np1)

  call radiation_boundary_condition(1, Nr, r, dr, dt, &
       scalar_nm1, Pi_nm1, scalar_n, Phi_n, a_n, alpha_n, &
       scalar_np1, Phi_np1, Pi_np1)

  Pi_np1(0) = -Pi_n(0) + Pi_np1(1) + Pi_n(1)
  call solve_metric(Nr, r, dr, Phi_np1, Pi_np1, a_np1, alpha_np1)
  t = t + 0.5_dp*dt
  call update_central_density_max(t, Phi_np1, Pi_np1, a_np1, rho_c_max, rho_c_max_time)

  scalar_n = scalar_np1
  Phi_n = Phi_np1
  Pi_n = Pi_np1
  a_n = a_np1
  alpha_n = alpha_np1

  do n = 2, Nt
     Phi_np1(0) = 0.0_dp
     a_np1(0) = 1.0_dp
     alpha_np1(0) = 1.0_dp

     call time_step(n, Nr, r, dr, dt, &
          scalar_n, Phi_n, Pi_n, a_n, alpha_n, &
          Phi_nm1, Pi_nm1, a_nm1, alpha_nm1, &
          Phi_np1, Pi_np1, scalar_np1)

     call radiation_boundary_condition(n, Nr, r, dr, dt, &
          scalar_nm1, Pi_nm1, scalar_n, Phi_n, a_n, alpha_n, &
          scalar_np1, Phi_np1, Pi_np1)

     Pi_np1(0) = -Pi_n(0) + Pi_np1(1) + Pi_n(1)
     call solve_metric(Nr, r, dr, Phi_np1, Pi_np1, a_np1, alpha_np1)
     t = t + dt
     call update_central_density_max(t, Phi_np1, Pi_np1, a_np1, rho_c_max, rho_c_max_time)

     if (has_nan(Nr, scalar_np1, Phi_np1, Pi_np1, a_np1, alpha_np1)) then
        error stop "La evolucion produjo NaN."
     end if

     if (minval(alpha_np1) < LAPSE_COLLAPSE_CRITERION) then
        lapse_collapsed = .true.
        write(*,'(A,I0,A,ES12.4)') "Colapso del lapse en n = ", n, ", t = ", t
        call output_snapshot(n, Nr, r, scalar_np1, Phi_np1, Pi_np1, a_np1, alpha_np1)
        write(*,'(A,I0,A,I0,A,ES12.4)') "Paso ", n, "/", Nt, ", t = ", t
        exit
     end if

     if (mod(n, OUTPUT_EVERY) == 0 .or. n == Nt) then
        call output_snapshot(n, Nr, r, scalar_np1, Phi_np1, Pi_np1, a_np1, alpha_np1)
        write(*,'(A,I0,A,I0,A,ES12.4)') "Paso ", n, "/", Nt, ", t = ", t
     end if

     scalar_nm1 = scalar_n
     scalar_n = scalar_np1
     Phi_nm1 = Phi_n
     Phi_n = Phi_np1
     Pi_nm1 = Pi_n
     Pi_n = Pi_np1
     a_nm1 = a_n
     a_n = a_np1
     alpha_nm1 = alpha_n
     alpha_n = alpha_np1
  end do

  if (lapse_collapsed) then
     write(*,*) "Resultado: Colapso de la funcion lapso."
  else
     write(*,*) "Resultado: evolucion terminada sin colapso de la funcion lapso."
  end if

  call output_central_density_max(phi0, phi0_input, phi0_decimals, rho_c_max)
  write(*,'(A)') "LEYENDA_RHO_MAX:"
  write(*,'(A,A)') "  phi0_input = ", trim(phi0_input)
  call print_decimal_value("  phi0_real = ", phi0, phi0_decimals)
  call print_decimal_value("  rho_c_max = ", rho_c_max, phi0_decimals)

contains

  subroutine delete_previous_outputs()

    integer :: exitstat
    integer :: cmdstat

    call execute_command_line( &
         "cmd /c if exist CEesferico1D_*.dat del /q CEesferico1D_*.dat", &
         wait=.true., exitstat=exitstat, cmdstat=cmdstat)

    if (cmdstat /= 0 .or. exitstat /= 0) then
       write(*,*) "Aviso: no se pudieron eliminar salidas anteriores CEesferico1D_*.dat."
    end if

    call execute_command_line( &         
         "cmd /c if exist mass_*.dat del /q mass_*.dat", &
         wait=.true., exitstat=exitstat, cmdstat=cmdstat)

    if (cmdstat /= 0 .or. exitstat /= 0) then
       write(*,*) "Aviso: no se pudieron eliminar salidas anteriores mass_*.dat."
    end if

    call execute_command_line( &
         "cmd /c if exist rho_central_max.dat del /q rho_central_max.dat", &
         wait=.true., exitstat=exitstat, cmdstat=cmdstat)

    if (cmdstat /= 0 .or. exitstat /= 0) then
       write(*,*) "Aviso: no se pudo eliminar salida anterior rho_central_max.dat."
    end if

  end subroutine delete_previous_outputs

  subroutine read_parameters(Nr, rmax, t_final, phi0, phi0_input, phi0_decimals)

    integer, intent(out) :: Nr
    real(dp), intent(out) :: rmax
    real(dp), intent(out) :: t_final
    real(dp), intent(out) :: phi0
    character(len=*), intent(out) :: phi0_input
    integer, intent(out) :: phi0_decimals
    character(len=100) :: arg

    if (command_argument_count() /= 4) then
       write(*,*) "Uso: CEesferico1D.exe Nr_intervalos rmax t_final phi0"
       write(*,*) "Ejemplo: CEesferico1D.exe 200 10.0 2.0 0.01"
       stop
    end if

    call get_command_argument(1, arg)
    read(arg,*) Nr
    Nr = Nr + 1
    call get_command_argument(2, arg)
    read(arg,*) rmax
    call get_command_argument(3, arg)
    read(arg,*) t_final
    call get_command_argument(4, arg)
    phi0_input = trim(arg)
    phi0_decimals = decimal_places(trim(arg))
    read(arg,*) phi0

    if (Nr < 4) error stop "Se requieren al menos tres intervalos radiales."
    if (rmax <= 0.0_dp) error stop "rmax debe ser positivo."
    if (t_final <= 0.0_dp) error stop "t_final debe ser positivo."

  end subroutine read_parameters

  integer function decimal_places(text)

    character(len=*), intent(in) :: text
    integer :: dot_pos
    integer :: exp_pos
    integer :: end_pos

    dot_pos = index(text, ".")
    if (dot_pos == 0) then
       decimal_places = 0
       return
    end if

    exp_pos = index(text, "E")
    if (exp_pos == 0) exp_pos = index(text, "e")

    if (exp_pos > 0) then
       end_pos = exp_pos - 1
    else
       end_pos = len_trim(text)
    end if

    decimal_places = max(0, end_pos - dot_pos)

  end function decimal_places

  subroutine print_decimal_value(label, value, decimals)

    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value
    integer, intent(in) :: decimals
    character(len=40) :: fmt

    write(fmt,'("(A,F0.",I0,")")') decimals
    write(*,fmt) label, value

  end subroutine print_decimal_value

  subroutine solve_metric(Nr, r, dr, Phi_field, Pi_field, a, alpha)

    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: Phi_field(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(inout) :: a(0:)
    real(dp), intent(inout) :: alpha(0:)
    integer :: j

    do j = 1, Nr - 1
       a(j) = Hamiltonian_constraint(j, r, dr, Phi_field, Pi_field, a)
       alpha(j) = polar_slicing_condition(j, r, dr, a, alpha)
    end do

    call rescaling_of_the_lapse(Nr, a, alpha)

  end subroutine solve_metric

  logical function has_nan(Nr, scalar, Phi_field, Pi_field, a, alpha)

    integer, intent(in) :: Nr
    real(dp), intent(in) :: scalar(0:)
    real(dp), intent(in) :: Phi_field(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(in) :: a(0:)
    real(dp), intent(in) :: alpha(0:)

    has_nan = any(ieee_is_nan(scalar(0:Nr-1))) .or. &
         any(ieee_is_nan(Phi_field(0:Nr-1))) .or. &
         any(ieee_is_nan(Pi_field(0:Nr-1))) .or. &
         any(ieee_is_nan(a(0:Nr-1))) .or. &
         any(ieee_is_nan(alpha(0:Nr-1)))

  end function has_nan

  real(dp) function energy_density(Phi_value, Pi_value, a_value)

    real(dp), intent(in) :: Phi_value
    real(dp), intent(in) :: Pi_value
    real(dp), intent(in) :: a_value

    energy_density = 0.5_dp*(Phi_value**2 + Pi_value**2)/a_value**2

  end function energy_density

  real(dp) function central_energy_density(Phi_field, Pi_field, a)

    real(dp), intent(in) :: Phi_field(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(in) :: a(0:)

    central_energy_density = energy_density(Phi_field(0), Pi_field(0), a(0))

  end function central_energy_density

  subroutine update_central_density_max(t, Phi_field, Pi_field, a, rho_c_max, rho_c_max_time)

    real(dp), intent(in) :: t
    real(dp), intent(in) :: Phi_field(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(in) :: a(0:)
    real(dp), intent(inout) :: rho_c_max
    real(dp), intent(inout) :: rho_c_max_time
    real(dp) :: rho_c

    rho_c = central_energy_density(Phi_field, Pi_field, a)
    if (rho_c > rho_c_max) then
       rho_c_max = rho_c
       rho_c_max_time = t
    end if

  end subroutine update_central_density_max

  real(dp) function misner_sharp_mass(r_value, a_value)

    real(dp), intent(in) :: r_value
    real(dp), intent(in) :: a_value

    misner_sharp_mass = 0.5_dp*r_value*(1.0_dp - 1.0_dp/a_value**2)

  end function misner_sharp_mass

  real(dp) function compactness(r_value, a_value)

    real(dp), intent(in) :: r_value
    real(dp), intent(in) :: a_value

    if (r_value > 0.0_dp) then
       compactness = 2.0_dp*misner_sharp_mass(r_value, a_value)/r_value
    else
       compactness = 0.0_dp
    end if

  end function compactness

  subroutine output_snapshot(n, Nr, r, scalar, Phi_field, Pi_field, a, alpha)

    integer, intent(in) :: n
    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: scalar(0:)
    real(dp), intent(in) :: Phi_field(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(in) :: a(0:)
    real(dp), intent(in) :: alpha(0:)
    integer :: j
    integer :: unit_number
    character(len=100) :: filename

    write(filename,'("CEesferico1D_",I8.8,".dat")') n
    open(newunit=unit_number, file=filename, status="replace", action="write")
    write(unit_number,'(A)') "# r scalar Phi Pi a alpha rho"
    do j = 0, Nr - 1
       write(unit_number,'(7ES24.15)') r(j), scalar(j), Phi_field(j), Pi_field(j), &
            a(j), alpha(j), energy_density(Phi_field(j), Pi_field(j), a(j))
    end do
    close(unit_number)

    call output_mass_aspect_function(n, Nr, r, a)

  end subroutine output_snapshot

  subroutine output_mass_aspect_function(n, Nr, r, a)

    integer, intent(in) :: n
    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: a(0:)
    integer :: j
    integer :: unit_number
    character(len=100) :: filename

    write(filename,'("mass_",I8.8,".dat")') n
    open(newunit=unit_number, file=filename, status="replace", action="write")
    write(unit_number,'(A)') "# r m_MS 2m_over_r"
    do j = 0, Nr - 1
       write(unit_number,'(3ES24.15)') r(j), misner_sharp_mass(r(j), a(j)), compactness(r(j), a(j))
    end do
    close(unit_number)

  end subroutine output_mass_aspect_function

  subroutine output_central_density_max(phi0, phi0_input, phi0_decimals, rho_c_max)

    real(dp), intent(in) :: phi0
    character(len=*), intent(in) :: phi0_input
    integer, intent(in) :: phi0_decimals
    real(dp), intent(in) :: rho_c_max
    integer :: unit_number
    character(len=40) :: fmt

    open(newunit=unit_number, file="rho_central_max.dat", status="replace", action="write")
    write(unit_number,'(A)') "# Leyenda de la corrida: densidad central maxima"
    write(unit_number,'(A)') "# rho = 0.5*(Phi**2 + Pi**2)/a**2 evaluada en r=0"
    write(unit_number,'(A,A)') "# phi0_input = ", trim(phi0_input)
    write(unit_number,'(A)') "# phi0_real rho_c_max"
    write(fmt,'("(F0.",I0,",1X,F0.",I0,")")') phi0_decimals, phi0_decimals
    write(unit_number,fmt) phi0, rho_c_max
    close(unit_number)

  end subroutine output_central_density_max

end program CEesferico1D
