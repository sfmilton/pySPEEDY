module gravity_wave_drag
    use types, only : p
    use params

    implicit none

    private
    public get_orographic_gravity_wave_drag

contains

    subroutine get_orographic_gravity_wave_drag(state, ug, vg, ut_gwd, vt_gwd)
        use model_state, only : ModelState_t
        use geometry, only : ModGeometry_t

        type(ModelState_t), intent(in), target :: state
        real(p), intent(in) :: ug(ix, il, kx), vg(ix, il, kx)
        real(p), intent(out) :: ut_gwd(ix, il, kx), vt_gwd(ix, il, kx)

        class(ModGeometry_t), pointer :: mod_geometry
        real(p) :: drag_tau, launch_sigma, oro_threshold, oro_scale
        real(p) :: mountain_factor, wind_factor, drag_coeff, profile
        real(p) :: low_level_wind
        integer :: i, j, k

        real(p), parameter :: min_tau_days = 1.0_p
        real(p), parameter :: wind_ref_speed = 20.0_p
        real(p), parameter :: max_wind_factor = 2.0_p

        mod_geometry => state%mod_geometry

        ut_gwd = 0.0_p
        vt_gwd = 0.0_p

        if (.not. state%orographic_gwd_enabled) return

        drag_tau = max(state%gwd_time_scale_days, min_tau_days) * 86400.0_p
        launch_sigma = min(max(state%gwd_launch_sigma, 0.05_p), 0.99_p)
        oro_threshold = max(state%gwd_oro_threshold_m, 0.0_p)
        oro_scale = max(state%gwd_oro_scale_m, 1.0_p)

        do i = 1, ix
            do j = 1, il
                mountain_factor = state%fmask_land(i, j) * &
                        min(max((state%orog(i, j) - oro_threshold) / oro_scale, 0.0_p), 1.0_p)
                if (mountain_factor <= 0.0_p) cycle

                low_level_wind = sqrt(ug(i, j, kx)**2 + vg(i, j, kx)**2)
                wind_factor = min(low_level_wind / wind_ref_speed, max_wind_factor)
                if (wind_factor <= 0.0_p) cycle

                drag_coeff = mountain_factor * wind_factor / drag_tau

                do k = 1, kx
                    if (mod_geometry%fsg(k) >= launch_sigma) cycle

                    ! Apply drag aloft with a smooth increase toward the model top.
                    profile = (launch_sigma - mod_geometry%fsg(k)) / launch_sigma
                    profile = profile * profile
                    ut_gwd(i, j, k) = -drag_coeff * profile * ug(i, j, k)
                    vt_gwd(i, j, k) = -drag_coeff * profile * vg(i, j, k)
                end do
            end do
        end do
    end subroutine

end module
