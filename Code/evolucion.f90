module evolution

  implicit none

  integer, parameter :: dp = kind(1.0d0)

contains

  !==========================================================================================
  ! CONDICION INICIAL GAUSSIANA EN COORDENADAS ESFERICAS en t=0 para phi, PHI, PI, a y alpha
  !==========================================================================================

  subroutine initial_condition(Nr, r, dr, phi0, r0, sigma, &
       scalar_old, Phi_field_old, Pi_field_old, a_old, alpha_old)

    implicit none

    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr, phi0, r0, sigma
    real(dp), intent(out) :: scalar_old(0:)
    real(dp), intent(out) :: Phi_field_old(0:)
    real(dp), intent(out) :: Pi_field_old(0:)
    real(dp), intent(out) :: a_old(0:)
    real(dp), intent(out) :: alpha_old(0:)

    integer :: j
    real(dp) :: factor, expfactor, gaussian

    do j = 0, Nr-1

       factor = (r(j) - r0)/(sigma**2)
       expfactor = (r(j) - r0)*factor
       gaussian = exp(-expfactor)

       ! phi(r,0) = phi0 exp(-(r-r0)^2/sigma^2)
       scalar_old(j) = phi0*gaussian

       ! Variable PI(r,0)= 0 Derivada temporal de la gaussiana.
       Pi_field_old(j) = 0.0_dp !sin velocidad inicial

      ! Variable PHI(r,0)= dphi / dr  Derivada radial de la gaussiana.

       if (j > 0) then

          Phi_field_old(j) = -2.0_dp*phi0*factor*gaussian

            ! Variable metrica a
          a_old(j) = Hamiltonian_constraint( &
               j, r, dr, Phi_field_old, Pi_field_old, a_old)

            ! Función lapso alpha
          alpha_old(j) = polar_slicing_condition( &
               j, r, dr, a_old, alpha_old)

       else

          ! Regularidad en el origen.
          Phi_field_old(j) = 0.0_dp
          a_old(j) = 1.0_dp
          alpha_old(j) = 1.0_dp

       end if

    end do

    call rescaling_of_the_lapse(Nr, a_old, alpha_old)

  end subroutine initial_condition

  !=========================================================
  ! EVOLUCION TEMPORAL DE phi, PHI Y PI
  !=========================================================

  subroutine time_step(n, Nr, r, dr, dt, &
       scalar_now, Phi_field_now, Pi_field_now, a_now, alpha_now, &
       Phi_field_old, Pi_field_old, a_old, alpha_old, &
       Phi_field_new, Pi_field_new, scalar_new)

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr, dt
    real(dp), intent(in) :: scalar_now(0:)
    real(dp), intent(in) :: Phi_field_now(0:)
    real(dp), intent(in) :: Pi_field_now(0:)
    real(dp), intent(in) :: a_now(0:)
    real(dp), intent(in) :: alpha_now(0:)
    real(dp), intent(in) :: Phi_field_old(0:)
    real(dp), intent(in) :: Pi_field_old(0:)
    real(dp), intent(in) :: a_old(0:)
    real(dp), intent(in) :: alpha_old(0:)
    real(dp), intent(out) :: Phi_field_new(0:)
    real(dp), intent(out) :: Pi_field_new(0:)
    real(dp), intent(inout) :: scalar_new(0:)

    integer :: j
    real(dp) :: inv_dr
    real(dp) :: Phi_Pi_dt_coef, phi_dt_coef
    real(dp) :: alpha_now_jm, alpha_old_j, alpha_now_j, alpha_now_jp
    real(dp) :: alpha_Pi_now_jm, alpha_Pi_old_j, alpha_Pi_now_j, alpha_Pi_now_jp
    real(dp) :: rhs_scalar, rhs_Phi, rhs_Pi
    real(dp) :: r_sqr_jm, r_sqr_jp, r_cbd_jm, r_cbd_jp
    real(dp) :: alpha_Phi_r2_over_a_jm, alpha_Phi_r2_over_a_jp

   !=========================================================
   ! COEFICIENTES TEMPORALES (separación entre los saltos)
   !=========================================================

    inv_dr = 1.0_dp/dr
   
   !Coeficientes temporales para Phi y Pi
    if (n >= 2) then
       Phi_Pi_dt_coef = 2.0_dp*dt
    else if (n == 1) then
       Phi_Pi_dt_coef = dt
    else
       Phi_Pi_dt_coef = 0.5_dp*dt
    end if

   !Coeficiente temporal para phi
    if (n >= 1) then
       phi_dt_coef = dt
    else
       phi_dt_coef = 0.5_dp*dt
    end if

   !Cálculos
    do j = 1, Nr-2
      
       alpha_now_jm = alpha_now(j-1)/a_now(j-1)
       alpha_old_j = alpha_old(j)/a_old(j)
       alpha_now_j = alpha_now(j)/a_now(j)
       alpha_now_jp = alpha_now(j+1)/a_now(j+1)

       alpha_Pi_now_jm = alpha_now_jm*Pi_field_now(j-1)
       alpha_Pi_old_j = alpha_old_j*Pi_field_old(j)
       alpha_Pi_now_j = alpha_now_j*Pi_field_now(j)
       alpha_Pi_now_jp = alpha_now_jp*Pi_field_now(j+1)


       !Definimos phi_new(0) para que evolucione a partir de j=1 con Adams-Boshforth
       if (j == 1) then
          rhs_scalar = 1.5_dp*alpha_now(0)*Pi_field_now(0)/a_now(0) - &
               0.5_dp*alpha_old(0)*Pi_field_old(0)/a_old(0)
          scalar_new(0) = scalar_now(0) + phi_dt_coef*rhs_scalar
       end if

       !A partir de j=2, evoluciona phi_new(j) con Adams-Boshforth

       rhs_scalar = 1.5_dp*alpha_Pi_now_j - 0.5_dp*alpha_Pi_old_j
       scalar_new(j) = scalar_now(j) + phi_dt_coef*rhs_scalar

       !Cálculos para PHI y PI con Adams-Boshforth
       r_sqr_jm = r(j-1)**2
       r_sqr_jp = r(j+1)**2
       r_cbd_jm = r_sqr_jm*r(j-1)
       r_cbd_jp = r_sqr_jp*r(j+1)
       alpha_Phi_r2_over_a_jm = r_sqr_jm*alpha_now_jm*Phi_field_now(j-1)
       alpha_Phi_r2_over_a_jp = r_sqr_jp*alpha_now_jp*Phi_field_now(j+1)

       !Calculo de PHI y PI nuevas con Adams-Boshforth

       rhs_Phi = 0.5_dp*inv_dr*(alpha_Pi_now_jp - alpha_Pi_now_jm)
       rhs_Pi = 3.0_dp*(alpha_Phi_r2_over_a_jp - alpha_Phi_r2_over_a_jm)/ (r_cbd_jp - r_cbd_jm)

       Phi_field_new(j) = Phi_field_old(j) + Phi_Pi_dt_coef*rhs_Phi
       Pi_field_new(j) = Pi_field_old(j) + Phi_Pi_dt_coef*rhs_Pi

    end do

  end subroutine time_step

  !=========================================================
  ! CONDICION DE FRONTERA DE RADIACION SALIENTE
  !=========================================================

  subroutine radiation_boundary_condition(n, Nr, r, dr, dt, scalar_old, Pi_field_old, &
       scalar_now, Phi_field_now, a_now, alpha_now, scalar_new, Phi_field_new, Pi_field_new)

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: Nr
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: scalar_old(0:)
    real(dp), intent(in) :: Pi_field_old(0:)
    real(dp), intent(in) :: scalar_now(0:)
    real(dp), intent(in) :: Phi_field_now(0:)
    real(dp), intent(in) :: a_now(0:)
    real(dp), intent(in) :: alpha_now(0:)
    real(dp), intent(inout) :: scalar_new(0:)
    real(dp), intent(inout) :: Phi_field_new(0:)
    real(dp), intent(inout) :: Pi_field_new(0:)

    integer :: J
    real(dp) :: inv_dr
    real(dp) :: tmp0
    real(dp) :: tmp1
    real(dp) :: tmp2
    real(dp) :: rhs_scalar
    real(dp) :: Phi_Pi_coef
    real(dp) :: rhs_Phi
    real(dp) :: r_sqd_Jm2
    real(dp) :: r_sqd_Jm1
    real(dp) :: r_sqd_J
    real(dp) :: r_cbd_Jm2
    real(dp) :: r_cbd_J
    real(dp) :: coef
    real(dp) :: term1
    real(dp) :: term2
    real(dp) :: term3
    real(dp) :: rhs_Pi

    inv_dr = 1.0_dp/dr
    J = Nr - 1

    tmp0 = -scalar_now(J)/r(J)
    tmp1 = 0.5_dp*inv_dr
    tmp2 = -tmp1*(3.0_dp*scalar_now(J) - 4.0_dp*scalar_now(J-1) + scalar_now(J-2))
    rhs_scalar = tmp0 + tmp2

    if (n >= 2) then
       Phi_Pi_coef = 2.0_dp*dt
    else if (n == 1) then
       Phi_Pi_coef = dt
    else
       Phi_Pi_coef = 0.5_dp*dt
    end if

    scalar_new(J) = scalar_old(J) + Phi_Pi_coef*rhs_scalar

    rhs_Phi = tmp1*(3.0_dp*scalar_new(J) - 4.0_dp*scalar_new(J-1) + scalar_new(J-2))
    Phi_field_new(J) = rhs_Phi

    r_sqd_Jm2 = r(J-2)**2
    r_sqd_Jm1 = r(J-1)**2
    r_sqd_J = r(J)**2
    r_cbd_Jm2 = r_sqd_Jm2*r(J-2)
    r_cbd_J = r_sqd_J*r(J)
    coef = 3.0_dp/(r_cbd_J - r_cbd_Jm2)

    term1 = 3.0_dp*r_sqd_J*(alpha_now(J)/a_now(J))*Phi_field_now(J)
    term2 = 4.0_dp*r_sqd_Jm1*(alpha_now(J-1)/a_now(J-1))*Phi_field_now(J-1)
    term3 = r_sqd_Jm2*(alpha_now(J-2)/a_now(J-2))*Phi_field_now(J-2)
    rhs_Pi = coef*(term1 - term2 + term3)

    Pi_field_new(J) = Pi_field_old(J) + Phi_Pi_coef*rhs_Pi

  end subroutine radiation_boundary_condition

  !=========================================================
  ! RESTRICCION HAMILTONIANA: CALCULO PUNTUAL DE a(j)
  !=========================================================

  function Hamiltonian_constraint( j, r, dr, Phi, Pi_field, a) result(valor_a)

    implicit none

    integer, intent(in) :: j
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: Phi(0:)
    real(dp), intent(in) :: Pi_field(0:)
    real(dp), intent(in) :: a(0:)
    real(dp) :: valor_a

    integer, parameter :: NEWTON_MAX_ITER = 300
    real(dp), parameter :: NEWTON_TOL = 1.0d-8
    real(dp), parameter :: NEWTON_MAX_DELTA = 0.25_dp
    real(dp) :: pi_const
    real(dp) :: inv_dr
    real(dp) :: A_prev
    real(dp) :: Phi_prom
    real(dp) :: Pi_prom
    real(dp) :: Phi_sqr
    real(dp) :: Pi_sqr
    real(dp) :: r_prom
    real(dp) :: Phi_Pi_term
    real(dp) :: inv_r_prom
    real(dp) :: A_old
    real(dp) :: A_new
    real(dp) :: tmp0
    real(dp) :: f
    real(dp) :: df
    real(dp) :: delta_A
    integer :: iter

    pi_const = acos(-1.0_dp)
    inv_dr = 1.0_dp/dr

    A_prev = log(a(j-1))
    Phi_prom = 0.5_dp*(Phi(j) + Phi(j-1)) ! Promedio de Phi en j y j-1 
    Pi_prom = 0.5_dp*(Pi_field(j) + Pi_field(j-1)) ! Promedio de Pi en j y j-1
    Phi_sqr = Phi_prom**2 ! Promedio de Phi al cuadrado en j y j-1
    Pi_sqr = Pi_prom**2 ! Promedio de Pi al cuadrado en j y j-1
    r_prom = 0.5_dp*(r(j) + r(j-1)) ! Promedio de r en j y j-1

    Phi_Pi_term = 2.0_dp*pi_const*r_prom*(Phi_sqr + Pi_sqr) !2*pi*r*(Phi^2 + Pi^2) en el punto medio entre j y j-1
    inv_r_prom = 0.5_dp/r_prom ! 1/2r

    A_old = A_prev
    A_new = A_old
    iter = 0

    do
       A_old = A_new

       tmp0 = inv_r_prom*exp(A_old + A_prev)
       f = inv_dr*(A_old - A_prev) + tmp0 - inv_r_prom - Phi_Pi_term
       df = inv_dr + tmp0

       delta_A = -f/df
       delta_A = max(-NEWTON_MAX_DELTA, min(NEWTON_MAX_DELTA, delta_A))
       A_new = A_old + delta_A   !Newton-Raphson amortiguado para encontrar f(A) = 0
       iter = iter + 1

       if (abs(A_new - A_old) <= NEWTON_TOL .or. iter > NEWTON_MAX_ITER) exit !Las iteraciones se detienen por tolerancia o por maximo de iteraciones.
    end do 

    if (iter > NEWTON_MAX_ITER) then
       write(*,*) "Newton warning: no convergio en j =", j, " iter =", iter
    end if

    valor_a = exp(A_new) !el código regresa: el valor de a(j) que satisface la restricción hamiltoniana en el punto Aj

  end function Hamiltonian_constraint

  !=========================================================
  ! POLAR SLICING: CALCULO PUNTUAL DE alpha(j)
  !=========================================================

  function polar_slicing_condition( j, r, dr, a, alpha) result(valor_alpha)

    implicit none

    integer, intent(in) :: j
    real(dp), intent(in) :: r(0:)
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: a(0:)
    real(dp), intent(in) :: alpha(0:)
    real(dp) :: valor_alpha

    real(dp) :: b
    real(dp) :: c
    real(dp) :: midway_r
    real(dp) :: d

    b = a(j) + a(j-1)
    c = a(j) - a(j-1)
    midway_r = 0.5_dp*(r(j) + r(j-1))

    d = (1.0_dp - 0.25_dp*b**2)/(2.0_dp*midway_r) - c/(dr*b)

    valor_alpha = alpha(j-1)*(1.0_dp - d*dr)/(1.0_dp + d*dr)

  end function polar_slicing_condition

  !=========================================================
  ! REESCALAMIENTO DEL LAPSO alpha (Velocidad máxima)
  !=========================================================

  subroutine rescaling_of_the_lapse(Nr, a, alpha)

    implicit none

    integer, intent(in) :: Nr
    real(dp), intent(in) :: a(0:)
    real(dp), intent(inout) :: alpha(0:)

    integer :: j
    real(dp) :: kappa
    real(dp) :: kappa_new

    kappa = a(0)/alpha(0)

    do j = 1, Nr-1
       kappa_new = a(j)/alpha(j)
       if (kappa_new < kappa) kappa = kappa_new
    end do

    do j = 0, Nr-1
       alpha(j) = alpha(j)*kappa
    end do

  end subroutine rescaling_of_the_lapse

end module evolution
