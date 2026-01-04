module rh_utils
  implicit none
  integer, parameter :: rk = selected_real_kind(12)

contains

  !--------------------------------------------------------
  ! Saturation vapor pressure (Pa)
  ! Unified Model / MATCH formulation
  ! T in Kelvin
  !--------------------------------------------------------
  real(rk) function esonT(T)
    implicit none
    real(rk), intent(in) :: T

    real(rk) :: log10e

    log10e = 10.79574_rk * (1.0_rk - 273.16_rk / T) &
           - 5.028_rk    * log10(T / 273.16_rk) &
           + 1.50475e-4_rk * (1.0_rk - 10.0_rk ** (-8.2969_rk * (T / 273.16_rk - 1.0_rk))) &
           + 0.42873e-3_rk * (10.0_rk ** (4.76955_rk * (1.0_rk - 273.16_rk / T)) - 1.0_rk) &
           + 0.78614_rk + 2.0_rk

    esonT = 10.0_rk ** log10e
  end function esonT

  !--------------------------------------------------------
  ! Relative humidity from mixing ratio, temperature, pressure
  !--------------------------------------------------------
  real(rk) function rh_vs_q_p_and_t(q, t_c, p)
    implicit none
    real(rk), intent(in) :: q     ! mixing ratio [kg/kg]
    real(rk), intent(in) :: t_c   ! temperature [C]
    real(rk), intent(in) :: p     ! pressure [Pa]

    real(rk) :: e, es

    e  = q * p / (0.622_rk + q)
    es = esonT(273.16_rk + t_c)

    rh_vs_q_p_and_t = 100.0_rk * e / es

    if (rh_vs_q_p_and_t > 100.0_rk) rh_vs_q_p_and_t = 100.0_rk
    if (rh_vs_q_p_and_t <   0.0_rk) rh_vs_q_p_and_t =   0.0_rk
  end function rh_vs_q_p_and_t

end module rh_utils
