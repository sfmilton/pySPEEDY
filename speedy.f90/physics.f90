module physics
    use types, only : p
    use params
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value

    implicit none

    private
    public get_physical_tendencies, reset_daily_tendency_diagnostics

contains

    !> Compute physical parametrization tendencies for u, v, t, q and add them
    !  to the dynamical grid-point tendencies
    subroutine get_physical_tendencies(state, j1, utend, vtend, ttend, qtend)
        use physical_constants, only : cp, grav, p0
        use sea_model, only : sea_coupling_flag
        use sppt, only : mu, gen_sppt
        use convection, only : get_convection_tendencies
        use large_scale_condensation, only : get_large_scale_condensation_tendencies
        use shortwave_radiation, only : get_shortwave_rad_fluxes, clouds
        use longwave_radiation, only : &
                get_downward_longwave_rad_fluxes, get_upward_longwave_rad_fluxes
        use surface_fluxes, only : get_surface_fluxes
        use vertical_diffusion, only : get_vertical_diffusion_tend
        use gravity_wave_drag, only : get_orographic_gravity_wave_drag
        use humidity, only : spec_hum_to_rel_hum
        use held_suarez, only : get_held_suarez_tendencies
        use physical_constants, only : alhc
        use model_state, only : ModelState_t
        use geometry, only : ModGeometry_t
        use Spectral, only : ModSpectral_t

        type(ModelState_t), intent(inout), target :: state
        integer, intent(in) :: j1

        real(p), intent(inout) :: utend(ix, il, kx) !! Zonal velocity tendency
        real(p), intent(inout) :: vtend(ix, il, kx) !! Meridional velocity tendency
        real(p), intent(inout) :: ttend(ix, il, kx) !! Temperature tendency
        real(p), intent(inout) :: qtend(ix, il, kx) !! Specific humidity tendency

        real(p), allocatable, dimension(:, :, :) :: tt_rlw
        real(p), allocatable, dimension(:, :, :) :: ut_pbl, vt_pbl, ut_gwd, vt_gwd
        real(p), allocatable, dimension(:, :, :) :: ut_hs, vt_hs
        real(p), allocatable, dimension(:, :, :) :: tt_pbl, qt_pbl
        real(p), allocatable, dimension(:, :, :) :: tt_sw, tt_lw, tt_hs

        integer :: i, j, k

        complex(p), allocatable, dimension(:, :) :: ucos, vcos
        real(p), allocatable, dimension(:, :) :: pslg, rps, gse
        real(p), allocatable, dimension(:, :) :: psg, ts, tskin, u0, v0, t0
        real(p), allocatable, dimension(:, :) :: cloudc, clstr, cltop, prtop, column_water_vapor
        logical, allocatable, dimension(:, :) :: has_cloud_top, has_conv_cloud_top
        integer, allocatable, dimension(:, :) :: iptop, icnv
        integer, allocatable :: icltop(:, :, :)

        real(p), allocatable, dimension(:, :, :) :: ug, vg, tg, qg, phig
        real(p), allocatable, dimension(:, :, :) :: utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
        real(p), allocatable, dimension(:, :, :) :: se, rh, qsat
        real(p), allocatable, dimension(:, :, :) :: tt_cnv, qt_cnv, tt_lsc, qt_lsc
        real(p), allocatable, dimension(:, :, :) :: sppt_pattern

        class(ModGeometry_t), pointer :: mod_geometry
        class(ModSpectral_t), pointer :: mod_spectral
        mod_geometry => state%mod_geometry
        mod_spectral => state%mod_spectral

        allocate(tt_cnv(ix, il, kx), qt_cnv(ix, il, kx), tt_lsc(ix, il, kx), qt_lsc(ix, il, kx))
        allocate(tt_rlw(ix, il, kx), ut_pbl(ix, il, kx), vt_pbl(ix, il, kx), tt_pbl(ix, il, kx))
        allocate(ut_gwd(ix, il, kx), vt_gwd(ix, il, kx), ut_hs(ix, il, kx), vt_hs(ix, il, kx))
        allocate(tt_sw(ix, il, kx), tt_lw(ix, il, kx), tt_hs(ix, il, kx))
        allocate(qt_pbl(ix, il, kx), sppt_pattern(ix, il, kx))
        allocate(ug(ix, il, kx), vg(ix, il, kx), tg(ix, il, kx), qg(ix, il, kx))
        allocate(phig(ix, il, kx), utend_dyn(ix, il, kx))
        allocate(vtend_dyn(ix, il, kx), ttend_dyn(ix, il, kx))
        allocate(qtend_dyn(ix, il, kx), se(ix, il, kx), rh(ix, il, kx), qsat(ix, il, kx))

        allocate(ucos(mx, nx), vcos(mx, nx))

        allocate(pslg(ix, il), rps(ix, il), gse(ix, il), psg(ix, il))
        allocate(ts(ix, il), tskin(ix, il), u0(ix, il), v0(ix, il))
        allocate(t0(ix, il), cloudc(ix, il), clstr(ix, il), cltop(ix, il), column_water_vapor(ix, il))
        allocate(prtop(ix, il), iptop(ix, il), icnv(ix, il))
        allocate(has_cloud_top(ix, il), has_conv_cloud_top(ix, il))

        allocate(icltop(ix, il, 2))

        ! Keep a copy of the original (dynamics only) tendencies
        utend_dyn = utend
        vtend_dyn = vtend
        ttend_dyn = ttend
        qtend_dyn = qtend
        tt_cnv = 0.0_p
        qt_cnv = 0.0_p
        tt_lsc = 0.0_p
        qt_lsc = 0.0_p
        tt_sw = 0.0_p
        tt_lw = 0.0_p
        ut_pbl = 0.0_p
        vt_pbl = 0.0_p
        ut_gwd = 0.0_p
        vt_gwd = 0.0_p
        ut_hs = 0.0_p
        vt_hs = 0.0_p
        tt_pbl = 0.0_p
        qt_pbl = 0.0_p
        tt_hs = 0.0_p
        cloudc = 0.0_p
        clstr = 0.0_p
        cltop = 0.0_p
        column_water_vapor = 0.0_p
        has_cloud_top = .false.
        has_conv_cloud_top = .false.

        ! =========================================================================
        ! Compute grid-point fields
        ! =========================================================================
        ! Convert model spectral variables to grid-point variables
        do k = 1, kx
            call mod_spectral%vort2vel(&
                    state%vor(:, :, k, j1), state%div(:, :, k, j1), ucos, vcos)

            ug(:, :, k) = mod_spectral%spec2grid(ucos, 2)
            vg(:, :, k) = mod_spectral%spec2grid(vcos, 2)
            tg(:, :, k) = mod_spectral%spec2grid(state%t(:, :, k, j1), 1)
            qg(:, :, k) = mod_spectral%spec2grid(state%tr(:, :, k, j1, 1), 1) ! q
            phig(:, :, k) = mod_spectral%spec2grid(state%phi(:, :, k), 1)

        end do

        pslg = mod_spectral%spec2grid(state%ps(:, :, j1), 1)

        ! =========================================================================
        ! Compute thermodynamic variables
        ! =========================================================================

        psg = exp(pslg)
        rps = 1.0 / psg

        if (state%held_suarez_mode) then
            call reset_dry_physics_outputs(state)
            call get_held_suarez_tendencies(&
                    psg, tg, ug, vg, mod_geometry%radang, mod_geometry%fsg, &
                    state%hs_trefc, state%hs_delta_ty, state%hs_delta_theta_z, state%hs_tmin, state%hs_sigma_b, &
                    state%hs_tau_a_days, state%hs_tau_s_days, state%hs_tau_f_days, &
                    state%hs_min_pressure_ratio, ut_hs, vt_hs, tt_hs)

            utend = utend + ut_hs
            vtend = vtend + vt_hs
            ttend = ttend + tt_hs
            qtend = 0.0
        else
            qg = max(qg, 0.0)
            se = cp * tg + phig
            column_water_vapor = 0.0_p
            do k = 1, kx
                column_water_vapor = column_water_vapor + qg(:, :, k) * 1.0e-3_p * p0 * psg * mod_geometry%dhs(k) / grav
            end do

            do k = 1, kx
                call spec_hum_to_rel_hum(tg(:, :, k), psg, mod_geometry%fsg(k), qg(:, :, k), &
                        rh(:, :, k), qsat(:, :, k))
            end do

            ! =========================================================================
            ! Precipitation
            ! =========================================================================

            ! Deep convection
            call get_convection_tendencies(psg, se, qg, qsat, iptop, state%cbmf, &
                    state%precnv, tt_cnv, qt_cnv, &
                    mod_geometry%fsg, mod_geometry%dhs, mod_geometry%wvi)

            do k = 2, kx
                tt_cnv(:, :, k) = tt_cnv(:, :, k) * rps * mod_geometry%grdscp(k)
                qt_cnv(:, :, k) = qt_cnv(:, :, k) * rps * mod_geometry%grdsig(k)
            end do

            icnv = kx - iptop

            ! Large-scale condensation
            call get_large_scale_condensation_tendencies(psg, qg, qsat, iptop, &
                    state%precls, tt_lsc, qt_lsc, mod_geometry%fsg, mod_geometry%dhs)

            ttend = ttend + tt_cnv + tt_lsc
            qtend = qtend + qt_cnv + qt_lsc

            ! =========================================================================
            ! Radiation (shortwave and longwave) and surface fluxes
            ! =========================================================================

            ! Since the shortwave tendencies may not computed at each time time state,
            ! the previous states are saved in the state%tt_rsw variable
            ! (Flux of short-wave radiation absorbed in each atmospheric layer).

            ! Compute shortwave tendencies and initialize lw transmissivity
            ! The shortwave radiation may be called at selected time steps
            if (state%compute_shortwave) then
                gse = (se(:, :, kx - 1) - se(:, :, kx)) / (phig(:, :, kx - 1) - phig(:, :, kx))

                call clouds(qg, rh, state%precnv, state%precls, iptop, gse, &
                        state%fmask_land, icltop, cloudc, clstr, state%qcloud_equiv)

                do i = 1, ix
                    do j = 1, il
                        cltop(i, j) = mod_geometry%sigh(icltop(i, j, 1) - 1) * p0 * psg(i, j) / 100.0_p
                        has_cloud_top(i, j) = cloudc(i, j) > 0.0_p
                        has_conv_cloud_top(i, j) = iptop(i, j) <= kx
                        if (has_conv_cloud_top(i, j)) then
                            prtop(i, j) = mod_geometry%sigh(iptop(i, j) - 1) * p0 * psg(i, j) / 100.0_p
                        else
                            prtop(i, j) = 0.0_p
                        end if
                    end do
                end do

                call get_shortwave_rad_fluxes(state, psg, qg, icltop, cloudc, clstr)

                do k = 1, kx
                    state%tt_rsw(:, :, k) = state%tt_rsw(:, :, k) * rps * mod_geometry%grdscp(k)
                end do
            end if
            tt_sw = state%tt_rsw

            ! Compute downward longwave fluxes
            call get_downward_longwave_rad_fluxes(&
                    tg, state%slrd, tt_rlw, state%fband, state%rad_flux, &
                    state%rad_tau2, state%rad_st4a, mod_geometry%wvi)

            ! Compute surface fluxes and land skin temperature
            call get_surface_fluxes(&
                    psg, ug, vg, tg, qg, rh, phig, &
                    state%phis0, state%fmask_land, state%forog, state%sst_am, &
                    & state%ssrd, state%slrd, state%ustr, state%vstr, &
                    state%shf, state%evap, state%slru, state%hfluxn, &
                    ts, tskin, u0, v0, t0, .true., &
                    state%alb_land, state%alb_sea, state%snowc, &
                    state%land_temp, state%soil_avail_water, &
                    mod_geometry%coa, mod_geometry%sigl, mod_geometry%wvi)

            ! Recompute sea fluxes in case of anomaly coupling
            if (sea_coupling_flag > 0) then
                call get_surface_fluxes(&
                        psg, ug, vg, tg, qg, rh, phig, state%phis0, state%fmask_land, state%forog, &
                        state%ssti_om, state%ssrd, state%slrd, &
                        state%ustr, state%vstr, state%shf, &
                        state%evap, state%slru, &
                        state%hfluxn, ts, tskin, u0, v0, t0, .false., &
                        state%alb_land, state%alb_sea, state%snowc, &
                        state%land_temp, state%soil_avail_water, &
                        mod_geometry%coa, mod_geometry%sigl, mod_geometry%wvi)
            end if

            ! Compute upward longwave fluxes, convert them to tendencies and add
            ! shortwave tendencies
            call get_upward_longwave_rad_fluxes(tg, ts, state%slrd, &
                    state%slru(:, :, 3), state%slr, &
                    state%olr, tt_rlw, state%fband, &
                    state%rad_flux, state%rad_tau2, state%rad_st4a, state%rad_strat_corr, &
                    mod_geometry%dhs)
            do k = 1, kx
                tt_rlw(:, :, k) = tt_rlw(:, :, k) * rps * mod_geometry%grdscp(k)
            end do
            tt_lw = tt_rlw

            ttend = ttend + tt_sw + tt_lw

            ! =========================================================================
            ! Planetary boundary later interactions with lower troposphere
            ! =========================================================================

            ! Vertical diffusion and shallow convection
            call get_vertical_diffusion_tend(&
                    se, rh, qg, qsat, phig, icnv, ut_pbl, vt_pbl, tt_pbl, qt_pbl, &
                    mod_geometry%fsg, mod_geometry%dhs, mod_geometry%sigh)

            ! Add tendencies due to surface fluxes
            ut_pbl(:, :, kx) = ut_pbl(:, :, kx) + state%ustr(:, :, 3) * rps * mod_geometry%grdsig(kx)
            vt_pbl(:, :, kx) = vt_pbl(:, :, kx) + state%vstr(:, :, 3) * rps * mod_geometry%grdsig(kx)
            tt_pbl(:, :, kx) = tt_pbl(:, :, kx) + state%shf(:, :, 3) * rps * mod_geometry%grdscp(kx)
            qt_pbl(:, :, kx) = qt_pbl(:, :, kx) + state%evap(:, :, 3) * rps * mod_geometry%grdsig(kx)

            utend = utend + ut_pbl
            vtend = vtend + vt_pbl
            ttend = ttend + tt_pbl
            qtend = qtend + qt_pbl

            call get_orographic_gravity_wave_drag(state, ug, vg, ut_gwd, vt_gwd)

            utend = utend + ut_gwd
            vtend = vtend + vt_gwd

            ! Add SPPT noise
            if (sppt_on) then
                sppt_pattern = gen_sppt(mod_spectral)

                ! The physical contribution to the tendency is *tend - *tend_dyn, where * is u, v, t, q
                do k = 1, kx
                    tt_cnv(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * tt_cnv(:, :, k)
                    tt_lsc(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * tt_lsc(:, :, k)
                    tt_sw(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * tt_sw(:, :, k)
                    tt_lw(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * tt_lw(:, :, k)
                    tt_pbl(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * tt_pbl(:, :, k)
                    qt_cnv(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * qt_cnv(:, :, k)
                    qt_lsc(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * qt_lsc(:, :, k)
                    qt_pbl(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * qt_pbl(:, :, k)
                    ut_pbl(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * ut_pbl(:, :, k)
                    vt_pbl(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * vt_pbl(:, :, k)
                    ut_gwd(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * ut_gwd(:, :, k)
                    vt_gwd(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * vt_gwd(:, :, k)
                    utend(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * (utend(:, :, k) - utend_dyn(:, :, k)) &
                            + utend_dyn(:, :, k)
                    vtend(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * (vtend(:, :, k) - vtend_dyn(:, :, k)) &
                            + vtend_dyn(:, :, k)
                    ttend(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * (ttend(:, :, k) - ttend_dyn(:, :, k)) &
                            + ttend_dyn(:, :, k)
                    qtend(:, :, k) = (1 + sppt_pattern(:, :, k) * mu(k)) * (qtend(:, :, k) - qtend_dyn(:, :, k)) &
                            + qtend_dyn(:, :, k)
                end do
            end if
        end if

        call update_daily_tendency_diagnostics( &
                state, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn, &
                utend - utend_dyn, vtend - vtend_dyn, ttend - ttend_dyn, qtend - qtend_dyn, &
                ut_pbl, vt_pbl, ut_gwd, vt_gwd, ut_hs, vt_hs, &
                tt_cnv, tt_lsc, tt_sw, tt_lw, tt_pbl, tt_hs, qt_cnv, qt_lsc, qt_pbl, &
                column_water_vapor, cloudc, clstr, cltop, prtop, has_cloud_top, has_conv_cloud_top, state%compute_shortwave)

        deallocate(ucos, vcos, pslg, rps, gse)
        deallocate(psg, ts, tskin, u0, v0, t0, cloudc, clstr, cltop, prtop, column_water_vapor)
        deallocate(has_cloud_top, has_conv_cloud_top)
        deallocate(iptop, icnv, icltop, ug, vg, tg, qg, phig, utend_dyn, vtend_dyn)
        deallocate(ttend_dyn, qtend_dyn, se, rh, qsat, sppt_pattern)
        deallocate(tt_cnv, qt_cnv, tt_lsc, qt_lsc)
        deallocate(tt_rlw, ut_pbl, vt_pbl, ut_gwd, vt_gwd, ut_hs, vt_hs, tt_pbl, qt_pbl)
        deallocate(tt_sw, tt_lw, tt_hs)
    end

    subroutine reset_dry_physics_outputs(state)
        use model_state, only : ModelState_t

        type(ModelState_t), intent(inout) :: state

        state%precnv = 0.0
        state%precls = 0.0
        state%snowcv = 0.0
        state%snowls = 0.0
        state%cbmf = 0.0
        state%tsr = 0.0
        state%ssrd = 0.0
        state%ssr = 0.0
        state%slrd = 0.0
        state%slr = 0.0
        state%olr = 0.0
        state%slru = 0.0
        state%ustr = 0.0
        state%vstr = 0.0
        state%shf = 0.0
        state%evap = 0.0
        state%hfluxn = 0.0
        state%tt_rsw = 0.0
        state%rad_flux = 0.0
        state%rad_tau2 = 0.0
        state%rad_st4a = 0.0
        state%rad_strat_corr = 0.0
        state%qcloud_equiv = 0.0
    end subroutine

    subroutine reset_daily_tendency_diagnostics(state)
        use model_state, only : ModelState_t

        type(ModelState_t), intent(inout) :: state

        state%ttend_dyn_accum = 0.0
        state%ttend_phy_accum = 0.0
        state%ttend_cnv_accum = 0.0
        state%ttend_lsc_accum = 0.0
        state%ttend_sw_accum = 0.0
        state%ttend_lw_accum = 0.0
        state%ttend_pbl_accum = 0.0
        state%ttend_hs_accum = 0.0
        state%qtend_dyn_accum = 0.0
        state%qtend_phy_accum = 0.0
        state%qtend_cnv_accum = 0.0
        state%qtend_lsc_accum = 0.0
        state%qtend_pbl_accum = 0.0
        state%utend_dyn_accum = 0.0
        state%vtend_dyn_accum = 0.0
        state%utend_phy_accum = 0.0
        state%vtend_phy_accum = 0.0
        state%utend_pbl_accum = 0.0
        state%vtend_pbl_accum = 0.0
        state%utend_gwd_accum = 0.0
        state%vtend_gwd_accum = 0.0
        state%utend_hs_accum = 0.0
        state%vtend_hs_accum = 0.0
        state%ustr_sfc_accum = 0.0
        state%vstr_sfc_accum = 0.0
        state%cloud_daily_count = 0
        state%cloud_cover_accum = 0.0
        state%stratiform_cloud_cover_accum = 0.0
        state%total_cloud_top_pressure_accum = 0.0
        state%total_cloud_top_count = 0.0
        state%conv_cloud_top_pressure_accum = 0.0
        state%conv_cloud_top_count = 0.0
        state%column_water_vapor_accum = 0.0
        state%precip_accum = 0.0
        state%evap_accum = 0.0
        state%toa_sw_down_accum = 0.0
        state%toa_sw_up_accum = 0.0
        state%toa_sw_net_accum = 0.0
        state%olr_accum = 0.0
        state%surface_lh_flux_accum = 0.0
        state%surface_sh_flux_accum = 0.0
        state%surface_sw_down_accum = 0.0
        state%surface_sw_up_accum = 0.0
        state%surface_sw_net_accum = 0.0
        state%surface_lw_down_accum = 0.0
        state%surface_lw_up_accum = 0.0
        state%surface_lw_net_accum = 0.0
        state%ttend_dyn_mean = 0.0
        state%ttend_phy_mean = 0.0
        state%ttend_cnv_mean = 0.0
        state%ttend_lsc_mean = 0.0
        state%ttend_sw_mean = 0.0
        state%ttend_lw_mean = 0.0
        state%ttend_pbl_mean = 0.0
        state%ttend_hs_mean = 0.0
        state%qtend_dyn_mean = 0.0
        state%qtend_phy_mean = 0.0
        state%qtend_cnv_mean = 0.0
        state%qtend_lsc_mean = 0.0
        state%qtend_pbl_mean = 0.0
        state%utend_dyn_mean = 0.0
        state%vtend_dyn_mean = 0.0
        state%utend_phy_mean = 0.0
        state%vtend_phy_mean = 0.0
        state%utend_pbl_mean = 0.0
        state%vtend_pbl_mean = 0.0
        state%utend_gwd_mean = 0.0
        state%vtend_gwd_mean = 0.0
        state%utend_hs_mean = 0.0
        state%vtend_hs_mean = 0.0
        state%ustr_sfc_mean = 0.0
        state%vstr_sfc_mean = 0.0
        state%cloud_cover_mean = 0.0
        state%stratiform_cloud_cover_mean = 0.0
        state%total_cloud_top_pressure_mean = 0.0
        state%conv_cloud_top_pressure_mean = 0.0
        state%column_water_vapor_mean = 0.0
        state%precip_mean = 0.0
        state%evap_mean = 0.0
        state%toa_sw_down_mean = 0.0
        state%toa_sw_up_mean = 0.0
        state%toa_sw_net_mean = 0.0
        state%olr_mean = 0.0
        state%surface_lh_flux_mean = 0.0
        state%surface_sh_flux_mean = 0.0
        state%surface_sw_down_mean = 0.0
        state%surface_sw_up_mean = 0.0
        state%surface_sw_net_mean = 0.0
        state%surface_lw_down_mean = 0.0
        state%surface_lw_up_mean = 0.0
        state%surface_lw_net_mean = 0.0
        state%ttend_daily_count = 0
    end subroutine

    subroutine update_daily_tendency_diagnostics( &
            state, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn, utend_phy, vtend_phy, ttend_phy, qtend_phy, &
            utend_pbl, vtend_pbl, utend_gwd, vtend_gwd, utend_hs, vtend_hs, &
            ttend_cnv, ttend_lsc, ttend_sw, ttend_lw, ttend_pbl, ttend_hs, qtend_cnv, qtend_lsc, qtend_pbl, &
            column_water_vapor, cloud_cover, stratiform_cloud_cover, total_cloud_top_pressure, conv_cloud_top_pressure, &
            has_cloud_top, has_conv_cloud_top, cloud_diag_available)
        use physical_constants, only : alhc
        use model_state, only : ModelState_t

        type(ModelState_t), intent(inout) :: state
        real(p), intent(in) :: utend_dyn(ix, il, kx), vtend_dyn(ix, il, kx)
        real(p), intent(in) :: ttend_dyn(ix, il, kx), qtend_dyn(ix, il, kx)
        real(p), intent(in) :: ttend_phy(ix, il, kx), qtend_phy(ix, il, kx)
        real(p), intent(in) :: utend_phy(ix, il, kx), vtend_phy(ix, il, kx)
        real(p), intent(in) :: utend_pbl(ix, il, kx), vtend_pbl(ix, il, kx)
        real(p), intent(in) :: utend_gwd(ix, il, kx), vtend_gwd(ix, il, kx)
        real(p), intent(in) :: utend_hs(ix, il, kx), vtend_hs(ix, il, kx)
        real(p), intent(in) :: ttend_cnv(ix, il, kx), ttend_lsc(ix, il, kx)
        real(p), intent(in) :: ttend_sw(ix, il, kx), ttend_lw(ix, il, kx)
        real(p), intent(in) :: ttend_pbl(ix, il, kx), ttend_hs(ix, il, kx)
        real(p), intent(in) :: qtend_cnv(ix, il, kx), qtend_lsc(ix, il, kx), qtend_pbl(ix, il, kx)
        real(p), intent(in) :: column_water_vapor(ix, il), cloud_cover(ix, il), stratiform_cloud_cover(ix, il)
        real(p), intent(in) :: total_cloud_top_pressure(ix, il), conv_cloud_top_pressure(ix, il)
        logical, intent(in) :: has_cloud_top(ix, il), has_conv_cloud_top(ix, il), cloud_diag_available
        real(p) :: day_scale, mean_scale, water_flux_scale

        if (mod(state%current_step, nsteps) == 0) then
            call reset_daily_tendency_diagnostics(state)
        end if

        state%utend_dyn_accum = state%utend_dyn_accum + utend_dyn
        state%vtend_dyn_accum = state%vtend_dyn_accum + vtend_dyn
        state%utend_phy_accum = state%utend_phy_accum + utend_phy
        state%vtend_phy_accum = state%vtend_phy_accum + vtend_phy
        state%utend_pbl_accum = state%utend_pbl_accum + utend_pbl
        state%vtend_pbl_accum = state%vtend_pbl_accum + vtend_pbl
        state%utend_gwd_accum = state%utend_gwd_accum + utend_gwd
        state%vtend_gwd_accum = state%vtend_gwd_accum + vtend_gwd
        state%utend_hs_accum = state%utend_hs_accum + utend_hs
        state%vtend_hs_accum = state%vtend_hs_accum + vtend_hs
        state%ustr_sfc_accum = state%ustr_sfc_accum + state%ustr(:, :, 3)
        state%vstr_sfc_accum = state%vstr_sfc_accum + state%vstr(:, :, 3)
        state%column_water_vapor_accum = state%column_water_vapor_accum + column_water_vapor
        if (cloud_diag_available) then
            state%cloud_daily_count = state%cloud_daily_count + 1
            state%cloud_cover_accum = state%cloud_cover_accum + cloud_cover
            state%stratiform_cloud_cover_accum = state%stratiform_cloud_cover_accum + stratiform_cloud_cover
            where (has_cloud_top)
                state%total_cloud_top_pressure_accum = state%total_cloud_top_pressure_accum + total_cloud_top_pressure
                state%total_cloud_top_count = state%total_cloud_top_count + 1.0_p
            end where
            where (has_conv_cloud_top)
                state%conv_cloud_top_pressure_accum = state%conv_cloud_top_pressure_accum + conv_cloud_top_pressure
                state%conv_cloud_top_count = state%conv_cloud_top_count + 1.0_p
            end where
        end if
        state%precip_accum = state%precip_accum + state%precnv + state%precls
        state%evap_accum = state%evap_accum + state%evap(:, :, 3)
        state%toa_sw_down_accum = state%toa_sw_down_accum + state%flux_solar_in
        state%toa_sw_up_accum = state%toa_sw_up_accum + (state%flux_solar_in - state%tsr)
        state%toa_sw_net_accum = state%toa_sw_net_accum + state%tsr
        state%olr_accum = state%olr_accum + state%olr
        state%surface_lh_flux_accum = state%surface_lh_flux_accum + alhc * state%evap(:, :, 3)
        state%surface_sh_flux_accum = state%surface_sh_flux_accum + state%shf(:, :, 3)
        state%surface_sw_down_accum = state%surface_sw_down_accum + state%ssrd
        state%surface_sw_up_accum = state%surface_sw_up_accum + (state%ssrd - state%ssr)
        state%surface_sw_net_accum = state%surface_sw_net_accum + state%ssr
        state%surface_lw_down_accum = state%surface_lw_down_accum + state%slrd
        state%surface_lw_up_accum = state%surface_lw_up_accum + state%slru(:, :, 3)
        state%surface_lw_net_accum = state%surface_lw_net_accum + (state%slrd - state%slru(:, :, 3))
        state%ttend_dyn_accum = state%ttend_dyn_accum + ttend_dyn
        state%ttend_phy_accum = state%ttend_phy_accum + ttend_phy
        state%ttend_cnv_accum = state%ttend_cnv_accum + ttend_cnv
        state%ttend_lsc_accum = state%ttend_lsc_accum + ttend_lsc
        state%ttend_sw_accum = state%ttend_sw_accum + ttend_sw
        state%ttend_lw_accum = state%ttend_lw_accum + ttend_lw
        state%ttend_pbl_accum = state%ttend_pbl_accum + ttend_pbl
        state%ttend_hs_accum = state%ttend_hs_accum + ttend_hs
        if (.not. state%held_suarez_mode) then
            state%qtend_dyn_accum = state%qtend_dyn_accum + qtend_dyn
            state%qtend_phy_accum = state%qtend_phy_accum + qtend_phy
            state%qtend_cnv_accum = state%qtend_cnv_accum + qtend_cnv
            state%qtend_lsc_accum = state%qtend_lsc_accum + qtend_lsc
            state%qtend_pbl_accum = state%qtend_pbl_accum + qtend_pbl
        end if
        state%ttend_daily_count = state%ttend_daily_count + 1

        day_scale = 86400.0_p / real(state%ttend_daily_count, p)
        mean_scale = 1.0_p / real(state%ttend_daily_count, p)
        water_flux_scale = 86.4_p / real(state%ttend_daily_count, p)
        state%utend_dyn_mean = state%utend_dyn_accum * day_scale
        state%vtend_dyn_mean = state%vtend_dyn_accum * day_scale
        state%utend_phy_mean = state%utend_phy_accum * day_scale
        state%vtend_phy_mean = state%vtend_phy_accum * day_scale
        state%utend_pbl_mean = state%utend_pbl_accum * day_scale
        state%vtend_pbl_mean = state%vtend_pbl_accum * day_scale
        state%utend_gwd_mean = state%utend_gwd_accum * day_scale
        state%vtend_gwd_mean = state%vtend_gwd_accum * day_scale
        state%utend_hs_mean = state%utend_hs_accum * day_scale
        state%vtend_hs_mean = state%vtend_hs_accum * day_scale
        state%ustr_sfc_mean = state%ustr_sfc_accum * mean_scale
        state%vstr_sfc_mean = state%vstr_sfc_accum * mean_scale
        state%column_water_vapor_mean = state%column_water_vapor_accum * mean_scale
        if (state%cloud_daily_count > 0) then
            state%cloud_cover_mean = state%cloud_cover_accum / real(state%cloud_daily_count, p)
            state%stratiform_cloud_cover_mean = state%stratiform_cloud_cover_accum / real(state%cloud_daily_count, p)
        else
            state%cloud_cover_mean = 0.0_p
            state%stratiform_cloud_cover_mean = 0.0_p
        end if
        where (state%total_cloud_top_count > 0.0_p)
            state%total_cloud_top_pressure_mean = &
                    state%total_cloud_top_pressure_accum / state%total_cloud_top_count
        elsewhere
            state%total_cloud_top_pressure_mean = ieee_value(state%total_cloud_top_pressure_mean, ieee_quiet_nan)
        end where
        where (state%conv_cloud_top_count > 0.0_p)
            state%conv_cloud_top_pressure_mean = &
                    state%conv_cloud_top_pressure_accum / state%conv_cloud_top_count
        elsewhere
            state%conv_cloud_top_pressure_mean = ieee_value(state%conv_cloud_top_pressure_mean, ieee_quiet_nan)
        end where
        state%precip_mean = state%precip_accum * water_flux_scale
        state%evap_mean = state%evap_accum * water_flux_scale
        state%toa_sw_down_mean = state%toa_sw_down_accum * mean_scale
        state%toa_sw_up_mean = state%toa_sw_up_accum * mean_scale
        state%toa_sw_net_mean = state%toa_sw_net_accum * mean_scale
        state%olr_mean = state%olr_accum * mean_scale
        state%surface_lh_flux_mean = state%surface_lh_flux_accum * mean_scale
        state%surface_sh_flux_mean = state%surface_sh_flux_accum * mean_scale
        state%surface_sw_down_mean = state%surface_sw_down_accum * mean_scale
        state%surface_sw_up_mean = state%surface_sw_up_accum * mean_scale
        state%surface_sw_net_mean = state%surface_sw_net_accum * mean_scale
        state%surface_lw_down_mean = state%surface_lw_down_accum * mean_scale
        state%surface_lw_up_mean = state%surface_lw_up_accum * mean_scale
        state%surface_lw_net_mean = state%surface_lw_net_accum * mean_scale
        state%ttend_dyn_mean = state%ttend_dyn_accum * day_scale
        state%ttend_phy_mean = state%ttend_phy_accum * day_scale
        state%ttend_cnv_mean = state%ttend_cnv_accum * day_scale
        state%ttend_lsc_mean = state%ttend_lsc_accum * day_scale
        state%ttend_sw_mean = state%ttend_sw_accum * day_scale
        state%ttend_lw_mean = state%ttend_lw_accum * day_scale
        state%ttend_pbl_mean = state%ttend_pbl_accum * day_scale
        state%ttend_hs_mean = state%ttend_hs_accum * day_scale
        state%qtend_dyn_mean = state%qtend_dyn_accum * day_scale
        state%qtend_phy_mean = state%qtend_phy_accum * day_scale
        state%qtend_cnv_mean = state%qtend_cnv_accum * day_scale
        state%qtend_lsc_mean = state%qtend_lsc_accum * day_scale
        state%qtend_pbl_mean = state%qtend_pbl_accum * day_scale
    end subroutine
end module
