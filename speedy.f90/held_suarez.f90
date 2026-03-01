!> author: Sean Milton (Codex assisted)
!  date: 01/03/2026
!> Held-Suarez dry benchmark forcing.
module held_suarez
    use types, only : p
    use params
    use physical_constants, only : akap

    implicit none

    private
    public get_held_suarez_tendencies

contains

    !> Compute Held-Suarez Newtonian cooling and Rayleigh drag tendencies.
    subroutine get_held_suarez_tendencies( &
            psa, tg, ug, vg, radang, sigma, trefc, delta_ty, delta_theta_z, tmin, sigma_b, &
            tau_a_days, tau_s_days, tau_f_days, min_pressure_ratio, utend, vtend, ttend)
        real(p), intent(in) :: psa(ix, il)
        real(p), intent(in) :: tg(ix, il, kx)
        real(p), intent(in) :: ug(ix, il, kx)
        real(p), intent(in) :: vg(ix, il, kx)
        real(p), intent(in) :: radang(il)
        real(p), intent(in) :: sigma(kx)
        real(p), intent(in) :: trefc, delta_ty, delta_theta_z, tmin, sigma_b
        real(p), intent(in) :: tau_a_days, tau_s_days, tau_f_days, min_pressure_ratio
        real(p), intent(out) :: utend(ix, il, kx)
        real(p), intent(out) :: vtend(ix, il, kx)
        real(p), intent(out) :: ttend(ix, il, kx)

        integer :: i, j, k
        real(p) :: sigma_weight, pressure_ratio, teq, sinphi, cosphi, kt, kv, k_a, k_s, k_f

        utend = 0.0_p
        vtend = 0.0_p
        ttend = 0.0_p
        k_a = 1.0_p / (tau_a_days * 86400.0_p)
        k_s = 1.0_p / (tau_s_days * 86400.0_p)
        k_f = 1.0_p / (tau_f_days * 86400.0_p)

        do k = 1, kx
            sigma_weight = max(0.0_p, (sigma(k) - sigma_b) / (1.0_p - sigma_b))

            do j = 1, il
                sinphi = sin(radang(j))
                cosphi = cos(radang(j))
                kt = k_a + (k_s - k_a) * sigma_weight * cosphi**4
                kv = k_f * sigma_weight

                do i = 1, ix
                    pressure_ratio = max(min_pressure_ratio, psa(i, j) * sigma(k))
                    teq = &
                        (trefc &
                        - delta_ty * sinphi**2 &
                        - delta_theta_z * log(pressure_ratio) * cosphi**2) &
                        * pressure_ratio**akap
                    teq = max(tmin, teq)

                    utend(i, j, k) = -kv * ug(i, j, k)
                    vtend(i, j, k) = -kv * vg(i, j, k)
                    ttend(i, j, k) = -kt * (tg(i, j, k) - teq)
                end do
            end do
        end do
    end subroutine

end module
