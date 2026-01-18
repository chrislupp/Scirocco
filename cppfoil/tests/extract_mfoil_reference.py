#!/usr/bin/env python3
"""
Extract reference data from mfoil for comprehensive C++ unit testing.
Saves data to HDF5 format for efficient loading in C++.
"""
import sys
sys.path.insert(0, '/Users/clupp/tools/analysis/ProjectTrysil/temp')
import mfoil
import numpy as np
import h5py

def extract_all_reference_data(filename):
    """Extract comprehensive reference data from mfoil and save to HDF5."""

    # Create NACA 0012 with ~130 panels
    M = mfoil.mfoil(naca='0012', npanel=130)
    M.setoper(alpha=0.0, Re=1e6)
    M.param.verb = 0

    with h5py.File(filename, 'w') as f:
        # ========== GEOMETRY DATA ==========
        geom = f.create_group('geometry')
        geom.create_dataset('N', data=M.foil.N)
        geom.create_dataset('x', data=M.foil.x)  # [2, N]
        geom.create_dataset('s', data=M.foil.s)  # [N]
        geom.create_dataset('t', data=M.foil.t)  # [2, N]

        # ========== PANEL INFLUENCE FUNCTIONS ==========
        panel = f.create_group('panel_influence')
        test_idx = 0
        for j in [10, 50, 100]:  # Panel indices
            Xj = M.foil.x[:, j:j+2]  # Panel endpoints
            for i in [5, 30, 70, 120]:  # Query points
                if i != j and i != j+1:
                    xi = M.foil.x[:, i]
                    a, b = mfoil.panel_linsource_stream(Xj, xi)
                    ax, bx = mfoil.panel_linsource_velocity(Xj, xi, np.array([1.0, 0.0]))
                    az, bz = mfoil.panel_linsource_velocity(Xj, xi, np.array([0.0, 1.0]))

                    grp = panel.create_group(f'test_{test_idx}')
                    grp.create_dataset('j', data=j)
                    grp.create_dataset('i', data=i)
                    grp.create_dataset('Xj', data=Xj)
                    grp.create_dataset('xi', data=xi)
                    grp.create_dataset('stream_a', data=a)
                    grp.create_dataset('stream_b', data=b)
                    grp.create_dataset('velocity_ax', data=ax)
                    grp.create_dataset('velocity_az', data=az)
                    grp.create_dataset('velocity_bx', data=bx)
                    grp.create_dataset('velocity_bz', data=bz)
                    test_idx += 1
        panel.attrs['num_tests'] = test_idx

        # ========== BUILD GAMMA (inviscid) ==========
        mfoil.build_gamma(M, M.oper.alpha)
        inviscid = f.create_group('inviscid')
        inviscid.create_dataset('gam', data=M.isol.gam)
        inviscid.create_dataset('gamref', data=M.isol.gamref)

        # ========== INIT THERMO ==========
        mfoil.init_thermo(M)

        # ========== WAKE ==========
        mfoil.build_wake(M)
        wake = f.create_group('wake')
        wake.create_dataset('N', data=M.wake.N)
        wake.create_dataset('x', data=M.wake.x)
        wake.create_dataset('s', data=M.wake.s)

        # ========== STAGNATION POINT ==========
        mfoil.stagpoint_find(M)
        stag = f.create_group('stagnation')
        stag.create_dataset('Istag', data=np.array(M.isol.Istag))
        stag.create_dataset('sstag', data=M.isol.sstag)
        stag.create_dataset('sstag_g', data=M.isol.sstag_g)

        # ========== XI and SURFACES ==========
        mfoil.identify_surfaces(M)
        mfoil.set_wake_gap(M)
        surfaces = f.create_group('surfaces')
        surfaces.create_dataset('Is_lower', data=np.array(list(M.vsol.Is[0])))
        surfaces.create_dataset('Is_upper', data=np.array(list(M.vsol.Is[1])))
        surfaces.create_dataset('Is_wake', data=np.array(list(M.vsol.Is[2])))
        surfaces.create_dataset('xi', data=M.isol.xi)
        surfaces.create_dataset('sgnue', data=M.isol.sgnue)
        surfaces.create_dataset('wgap', data=M.vsol.wgap)

        # ========== INVISCID UE ==========
        ue_inv = mfoil.get_ueinv(M)
        f.create_dataset('ue_inv', data=ue_inv)

        # ========== CALC UE_M ==========
        mfoil.calc_ue_m(M)
        ue_m = f.create_group('ue_m')
        ue_m.create_dataset('matrix', data=M.vsol.ue_m)
        # sigma_m and ue_sigma may be sparse, convert to dense
        if hasattr(M.vsol.sigma_m, 'toarray'):
            ue_m.create_dataset('sigma_m', data=M.vsol.sigma_m.toarray())
        else:
            ue_m.create_dataset('sigma_m', data=np.array(M.vsol.sigma_m))
        if hasattr(M.vsol.ue_sigma, 'toarray'):
            ue_m.create_dataset('ue_sigma', data=M.vsol.ue_sigma.toarray())
        else:
            ue_m.create_dataset('ue_sigma', data=np.array(M.vsol.ue_sigma))

        # ========== SAVE PARAM VALUES ==========
        # Save the actual param values used by mfoil
        param_grp = f.create_group('param')
        param_grp.create_dataset('Vinf', data=M.param.Vinf)
        param_grp.create_dataset('Minf', data=M.param.Minf)
        param_grp.create_dataset('muinf', data=M.param.muinf)
        param_grp.create_dataset('mu0', data=M.param.mu0)
        param_grp.create_dataset('rho0', data=M.param.rho0)
        param_grp.create_dataset('gam', data=M.param.gam)
        param_grp.create_dataset('KTb', data=M.param.KTb)
        param_grp.create_dataset('KTl', data=M.param.KTl)

        # ========== BL PROPERTY FUNCTIONS ==========
        bl_props = f.create_group('bl_properties')

        # Sample states: laminar, turbulent, near separation
        test_states = [
            {'th': 0.001, 'ds': 0.003, 'sa': 0.0, 'ue': 1.0, 'turb': False, 'wake': False, 'name': 'laminar1'},
            {'th': 0.002, 'ds': 0.004, 'sa': 0.0, 'ue': 0.8, 'turb': False, 'wake': False, 'name': 'laminar2'},
            {'th': 0.005, 'ds': 0.015, 'sa': 0.1, 'ue': 0.5, 'turb': True, 'wake': False, 'name': 'turb1'},
            {'th': 0.01, 'ds': 0.025, 'sa': 0.15, 'ue': 0.3, 'turb': True, 'wake': False, 'name': 'turb2'},
            {'th': 0.02, 'ds': 0.05, 'sa': 0.2, 'ue': 0.2, 'turb': True, 'wake': True, 'name': 'wake1'},
        ]

        for state in test_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param = mfoil.build_param(M, 0)
            param.turb = state['turb']
            param.wake = state['wake']

            grp = bl_props.create_group(state['name'])
            grp.create_dataset('U', data=U)
            grp.attrs['turb'] = state['turb']
            grp.attrs['wake'] = state['wake']

            # get_H
            H, H_U = mfoil.get_H(U)
            grp.create_dataset('H', data=H)
            grp.create_dataset('H_U', data=H_U)

            # get_Hk
            Hk, Hk_U = mfoil.get_Hk(U, param)
            grp.create_dataset('Hk', data=Hk)
            grp.create_dataset('Hk_U', data=Hk_U)

            # get_Hs
            Hs, Hs_U = mfoil.get_Hs(U, param)
            grp.create_dataset('Hs', data=Hs)
            grp.create_dataset('Hs_U', data=Hs_U)

            # get_Ret
            Ret, Ret_U = mfoil.get_Ret(U, param)
            grp.create_dataset('Ret', data=Ret)
            grp.create_dataset('Ret_U', data=Ret_U)

            # get_cf
            cf, cf_U = mfoil.get_cf(U, param)
            grp.create_dataset('cf', data=cf)
            grp.create_dataset('cf_U', data=cf_U)

            # get_cfxt at x=0.5
            x_test = 0.5
            cfxt, cfxt_U, cfxt_x = mfoil.get_cfxt(U, x_test, param)
            grp.create_dataset('cfxt', data=cfxt)
            grp.create_dataset('cfxt_U', data=cfxt_U)
            grp.create_dataset('cfxt_x', data=cfxt_x)

            # get_cDi
            cDi, cDi_U = mfoil.get_cDi(U, param)
            grp.create_dataset('cDi', data=cDi)
            grp.create_dataset('cDi_U', data=cDi_U)

            # get_cDixt at x=0.5
            cDixt, cDixt_U, cDixt_x = mfoil.get_cDixt(U, x_test, param)
            grp.create_dataset('cDixt', data=cDixt)
            grp.create_dataset('cDixt_U', data=cDixt_U)
            grp.create_dataset('cDixt_x', data=cDixt_x)

            # get_de
            de, de_U = mfoil.get_de(U, param)
            grp.create_dataset('de', data=de)
            grp.create_dataset('de_U', data=de_U)

            # get_Hss
            Hss, Hss_U = mfoil.get_Hss(U, param)
            grp.create_dataset('Hss', data=Hss)
            grp.create_dataset('Hss_U', data=Hss_U)

            # get_Hw (with wgap=0)
            Hw, Hw_U = mfoil.get_Hw(U, 0.0)
            grp.create_dataset('Hw', data=Hw)
            grp.create_dataset('Hw_U', data=Hw_U)

            if state['turb']:
                # get_Us
                Us, Us_U = mfoil.get_Us(U, param)
                grp.create_dataset('Us', data=Us)
                grp.create_dataset('Us_U', data=Us_U)

                # get_cteq
                cteq, cteq_U = mfoil.get_cteq(U, param)
                grp.create_dataset('cteq', data=cteq)
                grp.create_dataset('cteq_U', data=cteq_U)

                # get_cttr
                cttr, cttr_U = mfoil.get_cttr(U, param)
                grp.create_dataset('cttr', data=cttr)
                grp.create_dataset('cttr_U', data=cttr_U)

            if not state['turb']:
                # get_damp
                damp, damp_U = mfoil.get_damp(U, param)
                grp.create_dataset('damp', data=damp)
                grp.create_dataset('damp_U', data=damp_U)

        # ========== RESIDUAL STATION ==========
        residuals = f.create_group('residual_station')

        M.oper.viscous = True
        M.oper.initbl = True
        mfoil.init_boundary_layer(M)

        # Test stations on upper surface
        Is = list(M.vsol.Is[1])
        test_idx = 0
        for i in [2, 5, 10, 30]:
            if i < len(Is) - 1:
                idx1 = Is[i]
                idx2 = Is[i+1]

                xi_pair = np.array([M.isol.xi[idx1], M.isol.xi[idx2]])
                U_pair = M.glob.U[:, [idx1, idx2]]
                Aux = np.array([0.0, 0.0])

                param = mfoil.build_param(M, 1)
                param.turb = M.vsol.turb[idx2]
                param.simi = (i <= 1)

                R, R_U, R_x = mfoil.residual_station(param, xi_pair, U_pair, Aux)

                grp = residuals.create_group(f'test_{test_idx}')
                grp.create_dataset('surface', data=1)
                grp.create_dataset('i', data=i)
                grp.create_dataset('idx', data=np.array([idx1, idx2]))
                grp.create_dataset('xi', data=xi_pair)
                grp.create_dataset('U', data=U_pair)
                grp.create_dataset('Aux', data=Aux)
                grp.attrs['turb'] = bool(param.turb)
                grp.attrs['simi'] = bool(param.simi)
                grp.create_dataset('R', data=R)
                grp.create_dataset('R_U', data=R_U)
                grp.create_dataset('R_x', data=R_x)
                test_idx += 1
        residuals.attrs['num_tests'] = test_idx

        # ========== INITIAL BL STATE ==========
        bl_init = f.create_group('bl_init')
        bl_init.create_dataset('Nsys', data=M.glob.Nsys)
        bl_init.create_dataset('U', data=M.glob.U)
        bl_init.create_dataset('turb', data=np.array(M.vsol.turb, dtype=bool))

        # ========== BUILD GLOB SYS ==========
        mfoil.build_glob_sys(M)
        glob_sys = f.create_group('glob_sys')
        glob_sys.create_dataset('R', data=M.glob.R)
        glob_sys.create_dataset('R_norm', data=np.linalg.norm(M.glob.R))

        # ========== THWAITES INIT ==========
        thwaites_grp = f.create_group('thwaites')
        # Test various K and nu values
        K_vals = [10.0, 50.0, 100.0, 200.0]
        nu = M.param.muinf / M.param.rho0
        for i, K in enumerate(K_vals):
            th, ds = mfoil.thwaites_init(K, nu)
            grp = thwaites_grp.create_group(f'test_{i}')
            grp.create_dataset('K', data=K)
            grp.create_dataset('nu', data=nu)
            grp.create_dataset('theta', data=th)
            grp.create_dataset('delta_star', data=ds)
        thwaites_grp.attrs['num_tests'] = len(K_vals)

        # ========== STAGNATION STATE ==========
        stag_state = f.create_group('stagnation_state')
        # Test with first two points on upper surface
        Is_upper = list(M.vsol.Is[1])
        if len(Is_upper) >= 3:
            idx1 = Is_upper[1]
            idx2 = Is_upper[2]
            U_pair = M.glob.U[:, [idx1, idx2]]
            x_pair = np.array([M.isol.xi[idx1], M.isol.xi[idx2]])
            Ust, Ust_U, Ust_x, xst = mfoil.stagnation_state(U_pair, x_pair)
            stag_state.create_dataset('U', data=U_pair)
            stag_state.create_dataset('x', data=x_pair)
            stag_state.create_dataset('Ust', data=Ust)
            stag_state.create_dataset('Ust_U', data=Ust_U)
            stag_state.create_dataset('Ust_x', data=Ust_x)
            stag_state.create_dataset('xst', data=xst)

        # ========== WAKE SYS ==========
        wake_sys_grp = f.create_group('wake_sys')
        param = mfoil.build_param(M, 2)  # Wake surface
        R_wake, R_U_wake, J_wake = mfoil.wake_sys(M, param)
        wake_sys_grp.create_dataset('R', data=R_wake)
        wake_sys_grp.create_dataset('R_U', data=R_U_wake)
        wake_sys_grp.create_dataset('J', data=np.array(J_wake))

        # ========== RESIDUAL TRANSITION ==========
        res_tran = f.create_group('residual_transition')
        # Find the transition station on upper surface
        Is_upper = list(M.vsol.Is[1])
        tran_idx = None
        for i in range(1, len(Is_upper)):
            idx = Is_upper[i]
            idx_prev = Is_upper[i-1]
            xi_curr = M.isol.xi[idx]
            xi_prev = M.isol.xi[idx_prev]
            Xt = M.vsol.Xt[1, 0]  # Upper surface transition xi
            if xi_prev < Xt <= xi_curr:
                tran_idx = i
                break

        if tran_idx is not None:
            idx1 = Is_upper[tran_idx - 1]
            idx2 = Is_upper[tran_idx]
            xi_pair = np.array([M.isol.xi[idx1], M.isol.xi[idx2]])
            U_pair = M.glob.U[:, [idx1, idx2]]
            Aux = np.array([0.0, 0.0])
            param = mfoil.build_param(M, 1)
            R, R_U, R_x = mfoil.residual_transition(M, param, xi_pair, U_pair, Aux)
            res_tran.create_dataset('xi', data=xi_pair)
            res_tran.create_dataset('U', data=U_pair)
            res_tran.create_dataset('Aux', data=Aux)
            res_tran.create_dataset('Xt', data=M.vsol.Xt[1, 0])
            res_tran.create_dataset('R', data=R)
            res_tran.create_dataset('R_U', data=R_U)
            res_tran.create_dataset('R_x', data=R_x)
            res_tran.create_dataset('tran_idx', data=tran_idx)
            res_tran.create_dataset('idx1', data=idx1)
            res_tran.create_dataset('idx2', data=idx2)

        # ========== MARCH AMPLIFICATION ==========
        march_amp = f.create_group('march_amplification')
        # Store amplification factor (sa) for laminar points before transition
        Is_upper = list(M.vsol.Is[1])
        laminar_amp = []
        laminar_idx = []
        for i, idx in enumerate(Is_upper):
            if not M.vsol.turb[idx]:
                laminar_amp.append(M.glob.U[2, idx])
                laminar_idx.append(idx)
        march_amp.create_dataset('amplification', data=np.array(laminar_amp))
        march_amp.create_dataset('indices', data=np.array(laminar_idx))
        march_amp.create_dataset('ncrit', data=M.param.ncrit)

        # ========== NEWTON SOLVER ==========
        newton_grp = f.create_group('newton')

        # Store initial state before Newton iterations
        U_init = M.glob.U.copy()
        R_init = M.glob.R.copy()
        newton_grp.create_dataset('U_init', data=U_init)
        newton_grp.create_dataset('R_init', data=R_init)
        newton_grp.create_dataset('R_init_norm', data=np.linalg.norm(R_init))

        # Run a few Newton iterations and capture the data
        num_iter = min(5, M.param.niglob)
        for it in range(num_iter):
            iter_grp = newton_grp.create_group(f'iter_{it}')

            # Store state before solve
            iter_grp.create_dataset('U_before', data=M.glob.U.copy())

            # Solve for dU
            mfoil.solve_glob(M)
            iter_grp.create_dataset('dU', data=M.glob.dU.copy())

            # Update state (captures omega)
            # We need to capture omega, so let's manually compute it
            dU = M.glob.dU
            U = M.glob.U

            # Compute omega limits (matching mfoil logic)
            due_max = 0.25  # max relative ue change
            dth_max = 0.5   # max relative theta change
            dn_max = 2.0    # max amplification increase
            dct_max = 0.05  # max ctau change

            omega = 1.0
            Nsys = M.glob.Nsys

            # ue limit
            for i in range(Nsys):
                if abs(U[3,i]) > 1e-8:
                    om = due_max * abs(U[3,i]) / max(abs(dU[3,i]), 1e-12)
                    omega = min(omega, om)

            # theta limit
            for i in range(Nsys):
                if U[0,i] > 1e-10:
                    om = dth_max * U[0,i] / max(abs(dU[0,i]), 1e-12)
                    omega = min(omega, om)

            # amplification limit (laminar only)
            for i in range(Nsys):
                if not M.vsol.turb[i] and dU[2,i] > 0:
                    om = dn_max / max(dU[2,i], 1e-12)
                    omega = min(omega, om)

            # ctau limit (turbulent only)
            for i in range(Nsys):
                if M.vsol.turb[i] and U[2,i] > 1e-6:
                    om = dct_max * U[2,i] / max(abs(dU[2,i]), 1e-12)
                    omega = min(omega, om)

            iter_grp.create_dataset('omega', data=omega)

            # Actually update the state
            mfoil.update_state(M)
            iter_grp.create_dataset('U_after', data=M.glob.U.copy())

            # Update stagnation point
            mfoil.stagpoint_move(M)

            # Update transition
            mfoil.update_transition(M)

            # Rebuild global system
            mfoil.build_glob_sys(M)
            iter_grp.create_dataset('R_after', data=M.glob.R.copy())
            iter_grp.create_dataset('R_after_norm', data=np.linalg.norm(M.glob.R))

            # Check convergence
            if np.linalg.norm(M.glob.R) < M.param.rtol:
                newton_grp.attrs['converged_iter'] = it
                break

        newton_grp.attrs['num_iterations'] = min(it + 1, num_iter)

        # Store final converged state
        newton_grp.create_dataset('U_final', data=M.glob.U.copy())
        newton_grp.create_dataset('R_final', data=M.glob.R.copy())
        newton_grp.create_dataset('R_final_norm', data=np.linalg.norm(M.glob.R))

        # ========== SPLINE FUNCTIONS ==========
        spline_grp = f.create_group('spline')

        # QuadSeg - quadrature points and weights
        xq, wq = mfoil.quadseg()
        spline_grp.create_dataset('quadseg_x', data=xq)
        spline_grp.create_dataset('quadseg_w', data=wq)

        # Spline2D - create from foil points and evaluate
        # PP is {'X': CubicSpline, 'Y': CubicSpline} where CubicSpline has .x (breakpoints) and .c (coeffs)
        PP = mfoil.spline2d(M.foil.x)
        # Store spline coefficients - CubicSpline.x is knots, CubicSpline.c is coefficients
        spline_grp.create_dataset('PP_X_s', data=PP['X'].x)  # breakpoints (same as s values)
        spline_grp.create_dataset('PP_X_c', data=PP['X'].c)  # coefficients [4, N-1]
        spline_grp.create_dataset('PP_Y_s', data=PP['Y'].x)
        spline_grp.create_dataset('PP_Y_c', data=PP['Y'].c)

        # Evaluate spline at test points
        s_test = np.linspace(0, PP['X'].x[-1], 20)
        spline_grp.create_dataset('eval_s', data=s_test)
        xy_eval = mfoil.splineval(PP, s_test)
        spline_grp.create_dataset('eval_xy', data=xy_eval)

        # Evaluate tangent
        tan_eval = mfoil.splinetan(PP, s_test)
        spline_grp.create_dataset('eval_tan', data=tan_eval)

        # Evaluate derivatives (first and second)
        # Use scipy's CubicSpline.derivative() method
        DPX = PP['X'].derivative()
        DPY = PP['Y'].derivative()
        D2PX = DPX.derivative()  # Second derivative
        D2PY = DPY.derivative()

        deriv1_eval = np.vstack((DPX(s_test), DPY(s_test)))
        deriv2_eval = np.vstack((D2PX(s_test), D2PY(s_test)))
        spline_grp.create_dataset('eval_deriv1', data=deriv1_eval)
        spline_grp.create_dataset('eval_deriv2', data=deriv2_eval)

        # ========== PANEL INFLUENCE - VORTEX FUNCTIONS ==========
        panel_vortex = f.create_group('panel_vortex')
        test_idx = 0
        for j in [10, 50, 100]:  # Panel indices
            Xj = M.foil.x[:, j:j+2]  # Panel endpoints
            for i in [5, 30, 70, 120]:  # Query points
                if i != j and i != j+1:
                    xi = M.foil.x[:, i]

                    # Linear vortex stream
                    a_vs, b_vs = mfoil.panel_linvortex_stream(Xj, xi)

                    # Linear vortex velocity (with vdir)
                    a_vx, b_vx = mfoil.panel_linvortex_velocity(Xj, xi, np.array([1.0, 0.0]), False)
                    a_vz, b_vz = mfoil.panel_linvortex_velocity(Xj, xi, np.array([0.0, 1.0]), False)

                    # Const source stream
                    psi_cs = mfoil.panel_constsource_stream(Xj, xi)

                    # Const source velocity
                    v_cs_x = mfoil.panel_constsource_velocity(Xj, xi, np.array([1.0, 0.0]))
                    v_cs_z = mfoil.panel_constsource_velocity(Xj, xi, np.array([0.0, 1.0]))

                    grp = panel_vortex.create_group(f'test_{test_idx}')
                    grp.create_dataset('j', data=j)
                    grp.create_dataset('i', data=i)
                    grp.create_dataset('Xj', data=Xj)
                    grp.create_dataset('xi', data=xi)
                    # Vortex stream
                    grp.create_dataset('vortex_stream_a', data=a_vs)
                    grp.create_dataset('vortex_stream_b', data=b_vs)
                    # Vortex velocity
                    grp.create_dataset('vortex_vel_ax', data=a_vx)
                    grp.create_dataset('vortex_vel_bx', data=b_vx)
                    grp.create_dataset('vortex_vel_az', data=a_vz)
                    grp.create_dataset('vortex_vel_bz', data=b_vz)
                    # Const source stream
                    grp.create_dataset('const_source_stream', data=psi_cs)
                    # Const source velocity
                    grp.create_dataset('const_source_vel_x', data=v_cs_x)
                    grp.create_dataset('const_source_vel_z', data=v_cs_z)
                    test_idx += 1
        panel_vortex.attrs['num_tests'] = test_idx

        # ========== PANEL INFO ==========
        panel_info_grp = f.create_group('panel_info')
        test_idx = 0
        for j in [10, 50, 100]:
            Xj = M.foil.x[:, j:j+2]
            for i in [5, 30, 70, 120]:
                if i != j and i != j+1:
                    xi = M.foil.x[:, i]
                    # panel_info returns: t, n, x, z, d, r1, r2, theta1, theta2
                    t, n, x_loc, z_loc, d, r1, r2, theta1, theta2 = mfoil.panel_info(Xj, xi)

                    grp = panel_info_grp.create_group(f'test_{test_idx}')
                    grp.create_dataset('Xj', data=Xj)
                    grp.create_dataset('xi', data=xi)
                    grp.create_dataset('t', data=t)
                    grp.create_dataset('n', data=n)
                    grp.create_dataset('x_loc', data=x_loc)
                    grp.create_dataset('z_loc', data=z_loc)
                    grp.create_dataset('d', data=d)
                    grp.create_dataset('r1', data=r1)
                    grp.create_dataset('r2', data=r2)
                    grp.create_dataset('theta1', data=theta1)
                    grp.create_dataset('theta2', data=theta2)
                    test_idx += 1
        panel_info_grp.attrs['num_tests'] = test_idx

        # ========== INVISCID VELOCITY ==========
        inv_vel = f.create_group('inviscid_velocity')
        # Test at a few field points
        test_points = [
            np.array([0.5, 0.1]),   # Above airfoil
            np.array([0.5, -0.1]),  # Below airfoil
            np.array([1.5, 0.0]),   # In wake
            np.array([0.0, 0.2]),   # Near leading edge
        ]
        for idx, xi in enumerate(test_points):
            V, V_G = mfoil.inviscid_velocity(M.foil.x, M.isol.gam, M.param.Vinf,
                                             M.oper.alpha, xi, dolin=True)
            grp = inv_vel.create_group(f'test_{idx}')
            grp.create_dataset('xi', data=xi)
            grp.create_dataset('V', data=V)
            grp.create_dataset('V_G', data=V_G)
        inv_vel.attrs['num_tests'] = len(test_points)

        # ========== GET UE INV REF ==========
        ue_inv_ref = mfoil.get_ueinvref(M)
        f.create_dataset('ue_inv_ref', data=ue_inv_ref)

        # ========== ADDITIONAL BL PROPERTIES ==========
        bl_extra = f.create_group('bl_extra')

        # Test get_upw with two adjacent states
        for state_name in ['laminar1', 'turb1']:
            is_turb = 'turb' in state_name
            U1 = np.array([0.001, 0.003, 0.0 if not is_turb else 0.1, 1.0])
            U2 = np.array([0.0015, 0.004, 0.0 if not is_turb else 0.12, 0.9])

            param = mfoil.build_param(M, 0)
            param.turb = is_turb
            param.wake = False

            upw, upw_U = mfoil.get_upw(U1, U2, param)

            grp = bl_extra.create_group(f'upw_{state_name}')
            grp.create_dataset('U1', data=U1)
            grp.create_dataset('U2', data=U2)
            grp.create_dataset('upw', data=upw)
            grp.create_dataset('upw_U', data=upw_U)
            grp.attrs['turb'] = is_turb

        # Test get_cp
        test_ue = [0.5, 0.8, 1.0, 1.2]
        cp_grp = bl_extra.create_group('cp')
        param = mfoil.build_param(M, 0)
        for i, ue in enumerate(test_ue):
            cp, cp_U = mfoil.get_cp(ue, param)
            cp_grp.create_dataset(f'ue_{i}', data=ue)
            cp_grp.create_dataset(f'cp_{i}', data=cp)
            cp_grp.create_dataset(f'cp_U_{i}', data=cp_U)
        cp_grp.attrs['num_tests'] = len(test_ue)

        # Test get_Mach2
        mach2_grp = bl_extra.create_group('mach2')
        test_states_m2 = [
            {'U': np.array([0.001, 0.003, 0.0, 1.0]), 'name': 'm2_1'},
            {'U': np.array([0.002, 0.005, 0.1, 0.5]), 'name': 'm2_2'},
        ]
        for state in test_states_m2:
            m2, m2_U = mfoil.get_Mach2(state['U'], param)
            grp = mach2_grp.create_group(state['name'])
            grp.create_dataset('U', data=state['U'])
            grp.create_dataset('M2', data=m2)
            grp.create_dataset('M2_U', data=m2_U)

        # Test get_rho
        rho_grp = bl_extra.create_group('rho')
        for state in test_states_m2:
            rho, rho_U = mfoil.get_rho(state['U'], param)
            grp = rho_grp.create_group(state['name'])
            grp.create_dataset('U', data=state['U'])
            grp.create_dataset('rho', data=rho)
            grp.create_dataset('rho_U', data=rho_U)

        # Test get_uk (Karman-Tsien corrected speed)
        uk_grp = bl_extra.create_group('uk')
        test_ue_uk = [0.3, 0.5, 0.8, 1.0, 1.2]
        for i, ue in enumerate(test_ue_uk):
            uk, uk_U = mfoil.get_uk(ue, param)
            uk_grp.create_dataset(f'ue_{i}', data=ue)
            uk_grp.create_dataset(f'uk_{i}', data=uk)
            uk_grp.create_dataset(f'uk_U_{i}', data=uk_U)
        uk_grp.attrs['num_tests'] = len(test_ue_uk)

        # Test get_uq (shear lag equilibrium term)
        uq_grp = bl_extra.create_group('uq')
        # For get_uq we need intermediate values
        test_uq_states = [
            {'th': 0.005, 'ds': 0.015, 'sa': 0.1, 'ue': 0.5, 'turb': True, 'name': 'uq_turb1'},
            {'th': 0.01, 'ds': 0.025, 'sa': 0.15, 'ue': 0.3, 'turb': True, 'name': 'uq_turb2'},
        ]
        for state in test_uq_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param_uq = mfoil.build_param(M, 0)
            param_uq.turb = state['turb']
            param_uq.wake = False

            # Compute intermediate values
            ds, ds_U = state['ds'], np.array([0.0, 1.0, 0.0, 0.0])  # ds_U = [0, 1, 0, 0]
            cf, cf_U = mfoil.get_cf(U, param_uq)
            Hk, Hk_U = mfoil.get_Hk(U, param_uq)
            Ret, Ret_U = mfoil.get_Ret(U, param_uq)

            uq, uq_U = mfoil.get_uq(ds, ds_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param_uq)

            grp = uq_grp.create_group(state['name'])
            grp.create_dataset('U', data=U)
            grp.create_dataset('ds', data=ds)
            grp.create_dataset('ds_U', data=ds_U)
            grp.create_dataset('cf', data=cf)
            grp.create_dataset('cf_U', data=cf_U)
            grp.create_dataset('Hk', data=Hk)
            grp.create_dataset('Hk_U', data=Hk_U)
            grp.create_dataset('Ret', data=Ret)
            grp.create_dataset('Ret_U', data=Ret_U)
            grp.create_dataset('uq', data=uq)
            grp.create_dataset('uq_U', data=uq_U)
            grp.attrs['turb'] = state['turb']

        # Test get_cfutstag (cf*ue*theta at stagnation)
        cfut_grp = bl_extra.create_group('cfutstag')
        for state in test_uq_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param_stag = mfoil.build_param(M, 0)
            param_stag.turb = state['turb']
            param_stag.wake = False

            cfut, cfut_U = mfoil.get_cfutstag(U, param_stag)

            grp = cfut_grp.create_group(state['name'].replace('uq_', 'cfut_'))
            grp.create_dataset('U', data=U)
            grp.create_dataset('cfut', data=cfut)
            grp.create_dataset('cfut_U', data=cfut_U)
            grp.attrs['turb'] = state['turb']

        # Test get_cdutstag (cDi*ue*theta at stagnation)
        cdut_grp = bl_extra.create_group('cdutstag')
        for state in test_uq_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param_stag = mfoil.build_param(M, 0)
            param_stag.turb = state['turb']
            param_stag.wake = False

            cdut, cdut_U = mfoil.get_cdutstag(U, param_stag)

            grp = cdut_grp.create_group(state['name'].replace('uq_', 'cdut_'))
            grp.create_dataset('U', data=U)
            grp.create_dataset('cdut', data=cdut)
            grp.create_dataset('cdut_U', data=cdut_U)
            grp.attrs['turb'] = state['turb']

        # ========== TE INFO ==========
        te_info_grp = f.create_group('te_info')
        # TE_info returns: t, hTE, dtdx, tcp, tdp
        t_te, hTE, dtdx, tcp, tdp = mfoil.TE_info(M.foil.x)
        te_info_grp.create_dataset('t', data=t_te)
        te_info_grp.create_dataset('hTE', data=hTE)
        te_info_grp.create_dataset('dtdx', data=dtdx)
        te_info_grp.create_dataset('tcp', data=tcp)
        te_info_grp.create_dataset('tdp', data=tdp)

        # ========== INVISCID SOLVE + FORCES ==========
        # Create fresh mfoil for inviscid-only solve
        M_inv = mfoil.mfoil(naca='0012', npanel=130)
        M_inv.setoper(alpha=2.0, Re=0)  # Re=0 for inviscid
        M_inv.param.verb = 0
        mfoil.solve_inviscid(M_inv)  # Inviscid solve only
        mfoil.calc_force(M_inv)

        inv_forces = f.create_group('inviscid_forces')
        inv_forces.create_dataset('cl', data=M_inv.post.cl)
        inv_forces.create_dataset('cd', data=M_inv.post.cd)
        inv_forces.create_dataset('cm', data=M_inv.post.cm)
        inv_forces.create_dataset('cp', data=M_inv.post.cp)
        inv_forces.create_dataset('cpi', data=M_inv.post.cpi)
        inv_forces.create_dataset('alpha', data=M_inv.oper.alpha)

        # ========== FULL VISCOUS SOLVE + FORCES ==========
        # Create a fresh mfoil instance for full solve
        M2 = mfoil.mfoil(naca='0012', npanel=130)
        M2.setoper(alpha=0.0, Re=1e6)
        M2.param.verb = 0
        M2.solve()  # Full viscous solve

        forces_grp = f.create_group('forces')
        forces_grp.create_dataset('cl', data=M2.post.cl)
        forces_grp.create_dataset('cd', data=M2.post.cd)
        forces_grp.create_dataset('cm', data=M2.post.cm)
        forces_grp.create_dataset('cdpi', data=M2.post.cdpi)
        forces_grp.create_dataset('cdf', data=M2.post.cdf)

        # Store cp distribution
        forces_grp.create_dataset('cp', data=M2.post.cp)
        forces_grp.create_dataset('cpi', data=M2.post.cpi)

        # Store BL distributions
        dist_grp = forces_grp.create_group('distributions')
        dist_grp.create_dataset('th', data=M2.post.th)
        dist_grp.create_dataset('ds', data=M2.post.ds)
        dist_grp.create_dataset('sa', data=M2.post.sa)
        dist_grp.create_dataset('ue', data=M2.post.ue)
        dist_grp.create_dataset('cf', data=M2.post.cf)
        dist_grp.create_dataset('Ret', data=M2.post.Ret)
        dist_grp.create_dataset('Hk', data=M2.post.Hk)

        # Store converged state
        forces_grp.create_dataset('converged', data=M2.glob.conv)
        forces_grp.create_dataset('U_converged', data=M2.glob.U)

        # ========== BUILD_PARAM REFERENCE ==========
        build_param_grp = f.create_group('build_param')
        for si in [0, 1, 2]:
            param = mfoil.build_param(M, si)
            grp = build_param_grp.create_group(f'surface_{si}')
            grp.create_dataset('wake', data=param.wake)
            grp.create_dataset('turb', data=param.turb)
            grp.create_dataset('simi', data=param.simi)

        # ========== STATION_PARAM REFERENCE ==========
        station_param_grp = f.create_group('station_param')
        Is_lower = list(M.vsol.Is[0])
        Is_upper = list(M.vsol.Is[1])
        Istag = M.isol.Istag
        # Test indices: lower surface, upper surface laminar, upper turbulent, stagnation
        test_indices = []
        if len(Is_lower) > 2:
            test_indices.append(Is_lower[2])
        if len(Is_upper) > 5:
            test_indices.append(Is_upper[5])
        if len(Is_upper) > 0:
            test_indices.append(Is_upper[-1])
        test_indices.extend([Istag[0], Istag[1]])
        station_param_grp.create_dataset('test_indices', data=np.array(test_indices))
        for idx in test_indices:
            si = 0 if idx in Is_lower else 1
            param = mfoil.build_param(M, si)
            mfoil.station_param(M, param, idx)
            grp = station_param_grp.create_group(f'idx_{idx}')
            grp.create_dataset('turb', data=param.turb)
            grp.create_dataset('simi', data=param.simi)
            grp.create_dataset('wake', data=param.wake)

        # ========== MAKE_PANELS / SPLINE_CURVATURE REFERENCE ==========
        # Use M.geom.xpoint (raw NACA points) to test SplineCurvature
        panels_grp = f.create_group('make_panels')
        panels_grp.create_dataset('xpoint', data=M.geom.xpoint)
        # Panel output is already in geometry, but store explicitly
        panels_grp.create_dataset('x', data=M.foil.x)
        panels_grp.create_dataset('s', data=M.foil.s)
        panels_grp.create_dataset('t', data=M.foil.t)
        panels_grp.create_dataset('N', data=M.foil.N)

        # Direct spline_curvature call
        spline_curv_grp = f.create_group('spline_curvature')
        xpoint = M.geom.xpoint.copy()
        X, S, XS = mfoil.spline_curvature(xpoint, 131, 2.0, 0.1, None)  # None = automatic spacing
        spline_curv_grp.create_dataset('xpoint', data=xpoint)
        spline_curv_grp.create_dataset('N', data=131)
        spline_curv_grp.create_dataset('Ufac', data=2.0)
        spline_curv_grp.create_dataset('TEfac', data=0.1)
        spline_curv_grp.create_dataset('X', data=X)
        spline_curv_grp.create_dataset('S', data=S)
        spline_curv_grp.create_dataset('XS', data=XS)

        # ========== STAGNATION POINT MOVE REFERENCE ==========
        # Create fresh mfoil for stagnation point move test
        M_stag = mfoil.mfoil(naca='0012', npanel=130)
        M_stag.setoper(alpha=2.0, Re=1e6)  # Non-zero alpha to have stagnation point shift
        M_stag.param.verb = 0
        M_stag.oper.viscous = True
        M_stag.oper.initbl = True

        # Initialize through the standard sequence
        # Note: build_wake must come BEFORE stagpoint_find because
        # stagpoint_find extends xi to include wake points using wake.s
        mfoil.init_thermo(M_stag)
        mfoil.build_gamma(M_stag, M_stag.oper.alpha)
        mfoil.build_wake(M_stag)
        mfoil.stagpoint_find(M_stag)
        mfoil.identify_surfaces(M_stag)
        mfoil.set_wake_gap(M_stag)
        mfoil.calc_ue_m(M_stag)
        mfoil.init_boundary_layer(M_stag)
        mfoil.build_glob_sys(M_stag)

        stag_move = f.create_group('stag_point_move')
        stag_move.create_dataset('Istag_before', data=np.array(M_stag.isol.Istag))
        stag_move.create_dataset('sstag_before', data=M_stag.isol.sstag)
        stag_move.create_dataset('xi_before', data=M_stag.isol.xi.copy())

        # Do one Newton iteration to get a non-trivial state
        mfoil.solve_glob(M_stag)
        mfoil.update_state(M_stag)

        # Save state BEFORE stagpoint_move for exact reproduction
        stag_move.create_dataset('U_before_move', data=M_stag.glob.U.copy())
        stag_move.create_dataset('Istag_before_move', data=np.array(M_stag.isol.Istag))
        stag_move.create_dataset('sstag_before_move', data=M_stag.isol.sstag)

        # Now move stagnation point
        mfoil.stagpoint_move(M_stag)

        stag_move.create_dataset('Istag_after', data=np.array(M_stag.isol.Istag))
        stag_move.create_dataset('sstag_after', data=M_stag.isol.sstag)
        stag_move.create_dataset('xi_after', data=M_stag.isol.xi.copy())
        stag_move.create_dataset('alpha', data=M_stag.oper.alpha)

        # ========== UPDATE TRANSITION REFERENCE ==========
        # Use the same state we already have from stag_point_move
        update_tran = f.create_group('update_transition')
        update_tran.create_dataset('Xt_before', data=M_stag.vsol.Xt.copy())
        update_tran.create_dataset('turb_before', data=np.array(M_stag.vsol.turb, dtype=bool))
        # Export full state for exact reproduction
        update_tran.create_dataset('U_before', data=M_stag.glob.U.copy())
        update_tran.create_dataset('Istag', data=np.array(M_stag.isol.Istag))
        update_tran.create_dataset('sstag', data=M_stag.isol.sstag)
        update_tran.create_dataset('xi', data=M_stag.isol.xi.copy())
        update_tran.create_dataset('alpha', data=M_stag.oper.alpha)

        mfoil.update_transition(M_stag)

        update_tran.create_dataset('Xt_after', data=M_stag.vsol.Xt.copy())
        update_tran.create_dataset('turb_after', data=np.array(M_stag.vsol.turb, dtype=bool))
        update_tran.create_dataset('U_after', data=M_stag.glob.U.copy())

        # ========== GET DISTRIBUTIONS REFERENCE ==========
        # Use the fully converged solution from M2 (forces_grp already has cl, cd, etc.)
        dist_ref = f.create_group('get_distributions')
        dist_ref.create_dataset('x', data=M2.foil.x)
        dist_ref.create_dataset('th', data=M2.post.th)
        dist_ref.create_dataset('ds', data=M2.post.ds)
        dist_ref.create_dataset('sa', data=M2.post.sa)
        dist_ref.create_dataset('ue', data=M2.post.ue)
        dist_ref.create_dataset('cf', data=M2.post.cf)
        dist_ref.create_dataset('Ret', data=M2.post.Ret)
        dist_ref.create_dataset('Hk', data=M2.post.Hk)
        # Export full state for exact reproduction
        dist_ref.create_dataset('U', data=M2.glob.U.copy())
        dist_ref.create_dataset('Istag', data=np.array(M2.isol.Istag))
        dist_ref.create_dataset('sstag', data=M2.isol.sstag)
        dist_ref.create_dataset('xi', data=M2.isol.xi.copy())
        dist_ref.create_dataset('turb', data=np.array(M2.vsol.turb, dtype=bool))
        dist_ref.create_dataset('alpha', data=M2.oper.alpha)

        # ========== SOLVE VISCOUS / SOLVE COUPLED REFERENCE ==========
        # Create yet another fresh instance to capture iteration history
        M_visc = mfoil.mfoil(naca='0012', npanel=130)
        M_visc.setoper(alpha=0.0, Re=1e6)
        M_visc.param.verb = 0
        M_visc.oper.viscous = True
        M_visc.oper.initbl = True
        M_visc.param.niglob = 25
        M_visc.param.rtol = 1e-6

        solve_visc = f.create_group('solve_viscous')

        # Initialize (this is what solve_viscous does before solve_coupled)
        # Note: build_wake must come BEFORE stagpoint_find
        mfoil.init_thermo(M_visc)
        mfoil.build_gamma(M_visc, M_visc.oper.alpha)
        mfoil.build_wake(M_visc)
        mfoil.stagpoint_find(M_visc)
        mfoil.identify_surfaces(M_visc)
        mfoil.set_wake_gap(M_visc)
        mfoil.calc_ue_m(M_visc)
        mfoil.init_boundary_layer(M_visc)

        solve_visc.create_dataset('U_init', data=M_visc.glob.U.copy())
        solve_visc.create_dataset('ue_m', data=M_visc.vsol.ue_m.copy())

        # Run full coupled solve
        mfoil.solve_coupled(M_visc)

        solve_visc.create_dataset('U_final', data=M_visc.glob.U.copy())
        solve_visc.create_dataset('converged', data=M_visc.glob.conv)
        solve_visc.create_dataset('niter', data=M_visc.glob.niter if hasattr(M_visc.glob, 'niter') else 0)

        # Calculate and store forces
        mfoil.calc_force(M_visc)
        solve_visc.create_dataset('cl', data=M_visc.post.cl)
        solve_visc.create_dataset('cd', data=M_visc.post.cd)
        solve_visc.create_dataset('cm', data=M_visc.post.cm)

        # ========== GEOMETRY MODIFICATION FUNCTIONS REFERENCE ==========
        # Test mgeom_flap
        M_flap = mfoil.mfoil(naca='0012', npanel=130)
        M_flap.param.verb = 0
        geom_flap = f.create_group('geom_flap')
        geom_flap.create_dataset('xpoint_before', data=M_flap.geom.xpoint.copy())

        xzhinge = np.array([0.75, 0.0])  # Hinge at 75% chord
        eta = 10.0  # 10 degrees down
        mfoil.mgeom_flap(M_flap, xzhinge, eta)

        geom_flap.create_dataset('xpoint_after', data=M_flap.geom.xpoint.copy())
        geom_flap.create_dataset('x_after', data=M_flap.foil.x.copy())
        geom_flap.create_dataset('xzhinge', data=xzhinge)
        geom_flap.create_dataset('eta', data=eta)

        # Test mgeom_addcamber
        M_camb = mfoil.mfoil(naca='0012', npanel=130)
        M_camb.param.verb = 0
        geom_camb = f.create_group('geom_addcamber')
        geom_camb.create_dataset('xpoint_before', data=M_camb.geom.xpoint.copy())

        # Create a simple parabolic camber increment
        x_camb = np.linspace(0, 1, 21)
        z_camb = 0.02 * x_camb * (1 - x_camb)  # Max 0.5% camber at mid-chord
        xzcamb = np.vstack([x_camb, z_camb])
        mfoil.mgeom_addcamber(M_camb, xzcamb)

        geom_camb.create_dataset('xpoint_after', data=M_camb.geom.xpoint.copy())
        geom_camb.create_dataset('x_after', data=M_camb.foil.x.copy())
        geom_camb.create_dataset('xzcamb', data=xzcamb)

        # Test mgeom_derotate
        # First create an airfoil that is rotated (use cambered airfoil)
        M_derot = mfoil.mfoil(naca='2412', npanel=130)
        M_derot.param.verb = 0
        geom_derot = f.create_group('geom_derotate')
        geom_derot.create_dataset('xpoint_before', data=M_derot.geom.xpoint.copy())

        mfoil.mgeom_derotate(M_derot)

        geom_derot.create_dataset('xpoint_after', data=M_derot.geom.xpoint.copy())
        geom_derot.create_dataset('x_after', data=M_derot.foil.x.copy())

        # ========== REBUILD ISOL REFERENCE ==========
        # Test that gamma can be rebuilt from gamref at different alpha
        rebuild_isol = f.create_group('rebuild_isol')

        # Create fresh mfoil and build gamma at alpha=0
        M_reb = mfoil.mfoil(naca='0012', npanel=130)
        M_reb.setoper(alpha=0.0, Re=1e6)
        M_reb.param.verb = 0
        mfoil.build_gamma(M_reb, 0.0)  # Build at alpha=0 to get gamref

        rebuild_isol.create_dataset('gamref', data=M_reb.isol.gamref.copy())
        rebuild_isol.create_dataset('gam_alpha0', data=M_reb.isol.gam.copy())

        # Now test rebuilding at alpha=5
        test_alpha = 5.0
        M_reb.oper.alpha = test_alpha
        cosa = np.cos(np.radians(test_alpha))
        sina = np.sin(np.radians(test_alpha))
        gam_rebuilt = M_reb.isol.gamref[:, 0] * cosa + M_reb.isol.gamref[:, 1] * sina

        rebuild_isol.create_dataset('test_alpha', data=test_alpha)
        rebuild_isol.create_dataset('gam_alpha5', data=gam_rebuilt)

        # Also test at alpha=-3
        test_alpha2 = -3.0
        cosa2 = np.cos(np.radians(test_alpha2))
        sina2 = np.sin(np.radians(test_alpha2))
        gam_rebuilt2 = M_reb.isol.gamref[:, 0] * cosa2 + M_reb.isol.gamref[:, 1] * sina2

        rebuild_isol.create_dataset('test_alpha2', data=test_alpha2)
        rebuild_isol.create_dataset('gam_alpha_neg3', data=gam_rebuilt2)

        # ========== SOLVE INVISCID DEDICATED TEST ==========
        # Create fresh mfoil for inviscid solve reference data
        solve_inv = f.create_group('solve_inviscid')
        M_sinv = mfoil.mfoil(naca='0012', npanel=130)
        M_sinv.setoper(alpha=2.0, Re=0)  # Re=0 for inviscid
        M_sinv.param.verb = 0
        mfoil.solve_inviscid(M_sinv)

        solve_inv.create_dataset('gamma', data=M_sinv.isol.gam)
        solve_inv.create_dataset('gamref', data=M_sinv.isol.gamref)
        solve_inv.create_dataset('sgnue', data=M_sinv.isol.sgnue)
        solve_inv.create_dataset('alpha', data=M_sinv.oper.alpha)
        # AIC matrix is built internally during solve - we can verify by rebuilding gamma
        cosa = np.cos(np.radians(2.0))
        sina = np.sin(np.radians(2.0))
        gam_check = M_sinv.isol.gamref[:, 0] * cosa + M_sinv.isol.gamref[:, 1] * sina
        solve_inv.create_dataset('gam_from_gamref', data=gam_check)

        # ========== GET CFXT DEDICATED TEST ==========
        # Test GetCfxt with explicit laminar and turbulent states
        cfxt_grp = f.create_group('get_cfxt')

        # Laminar test cases
        lam_states = [
            {'th': 0.001, 'ds': 0.003, 'sa': 0.0, 'ue': 1.0, 'x': 0.3, 'name': 'laminar_1'},
            {'th': 0.0015, 'ds': 0.004, 'sa': 0.5, 'ue': 0.9, 'x': 0.5, 'name': 'laminar_2'},
        ]
        for state in lam_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param = mfoil.build_param(M, 0)
            param.turb = False
            param.wake = False
            cfxt, cfxt_U, cfxt_x = mfoil.get_cfxt(U, state['x'], param)
            grp = cfxt_grp.create_group(state['name'])
            grp.create_dataset('U', data=U)
            grp.create_dataset('x', data=state['x'])
            grp.create_dataset('cfxt', data=cfxt)
            grp.create_dataset('cfxt_U', data=cfxt_U)
            grp.create_dataset('cfxt_x', data=cfxt_x)
            grp.attrs['turb'] = False

        # Turbulent test cases
        turb_states = [
            {'th': 0.005, 'ds': 0.015, 'sa': 0.1, 'ue': 0.5, 'x': 0.5, 'name': 'turb_1'},
            {'th': 0.01, 'ds': 0.025, 'sa': 0.15, 'ue': 0.3, 'x': 0.7, 'name': 'turb_2'},
        ]
        for state in turb_states:
            U = np.array([state['th'], state['ds'], state['sa'], state['ue']])
            param = mfoil.build_param(M, 0)
            param.turb = True
            param.wake = False
            cfxt, cfxt_U, cfxt_x = mfoil.get_cfxt(U, state['x'], param)
            grp = cfxt_grp.create_group(state['name'])
            grp.create_dataset('U', data=U)
            grp.create_dataset('x', data=state['x'])
            grp.create_dataset('cfxt', data=cfxt)
            grp.create_dataset('cfxt_U', data=cfxt_U)
            grp.create_dataset('cfxt_x', data=cfxt_x)
            grp.attrs['turb'] = True

        # ========== UPWIND DEDICATED TEST ==========
        upwind_grp = f.create_group('upwind')
        # Test Upwind function with explicit states
        upwind_tests = [
            {'U1': [0.001, 0.003, 0.0, 1.0], 'U2': [0.0012, 0.0035, 0.0, 0.95], 'turb': False, 'name': 'laminar'},
            {'U1': [0.005, 0.015, 0.1, 0.5], 'U2': [0.006, 0.018, 0.12, 0.45], 'turb': True, 'name': 'turbulent'},
        ]
        for test in upwind_tests:
            U1 = np.array(test['U1'])
            U2 = np.array(test['U2'])
            param = mfoil.build_param(M, 0)
            param.turb = test['turb']
            param.wake = False
            upw, upw_U = mfoil.get_upw(U1, U2, param)
            grp = upwind_grp.create_group(test['name'])
            grp.create_dataset('U1', data=U1)
            grp.create_dataset('U2', data=U2)
            grp.create_dataset('upw', data=upw)
            grp.create_dataset('upw_U', data=upw_U)
            grp.attrs['turb'] = test['turb']

        # ========== EXTEND XI TO WAKE DEDICATED TEST ==========
        xi_wake_grp = f.create_group('extend_xi_to_wake')
        # Create fresh mfoil for this test
        M_xi = mfoil.mfoil(naca='0012', npanel=130)
        M_xi.setoper(alpha=0.0, Re=1e6)
        M_xi.param.verb = 0
        mfoil.init_thermo(M_xi)
        mfoil.build_gamma(M_xi, M_xi.oper.alpha)
        mfoil.build_wake(M_xi)
        mfoil.stagpoint_find(M_xi)
        # Save state BEFORE extending xi
        xi_before = M_xi.foil.s.copy()  # s values before extension
        N_before = M_xi.foil.N
        # Extend xi
        # In mfoil, xi is extended implicitly - we need to call identify_surfaces which does it
        mfoil.identify_surfaces(M_xi)
        xi_wake_grp.create_dataset('xi', data=M_xi.isol.xi)
        xi_wake_grp.create_dataset('N_foil', data=M_xi.foil.N)
        xi_wake_grp.create_dataset('N_wake', data=M_xi.wake.N)
        xi_wake_grp.create_dataset('foil_s', data=M_xi.foil.s)
        xi_wake_grp.create_dataset('wake_s', data=M_xi.wake.s)
        xi_wake_grp.create_dataset('sstag', data=M_xi.isol.sstag)

        # ========== IDENTIFY SURFACES DEDICATED TEST ==========
        surfaces_grp = f.create_group('identify_surfaces')
        surfaces_grp.create_dataset('Is_lower', data=np.array(list(M_xi.vsol.Is[0])))
        surfaces_grp.create_dataset('Is_upper', data=np.array(list(M_xi.vsol.Is[1])))
        surfaces_grp.create_dataset('Is_wake', data=np.array(list(M_xi.vsol.Is[2])))
        surfaces_grp.create_dataset('Istag', data=np.array(M_xi.isol.Istag))

        # ========== SET WAKE GAP DEDICATED TEST ==========
        wgap_grp = f.create_group('set_wake_gap')
        mfoil.set_wake_gap(M_xi)
        wgap_grp.create_dataset('wgap', data=M_xi.vsol.wgap)
        wgap_grp.create_dataset('Is_wake', data=np.array(list(M_xi.vsol.Is[2])))

        # ========== STORE TRANSITION DEDICATED TEST ==========
        store_tran = f.create_group('store_transition')
        # Initialize BL to get transition data
        mfoil.calc_ue_m(M_xi)
        M_xi.oper.viscous = True
        M_xi.oper.initbl = True
        mfoil.init_boundary_layer(M_xi)
        # Get transition location after BL initialization
        store_tran.create_dataset('Xt', data=M_xi.vsol.Xt)
        store_tran.create_dataset('turb', data=np.array(M_xi.vsol.turb, dtype=bool))
        store_tran.create_dataset('Is_lower', data=np.array(list(M_xi.vsol.Is[0])))
        store_tran.create_dataset('Is_upper', data=np.array(list(M_xi.vsol.Is[1])))
        store_tran.create_dataset('ncrit', data=M_xi.param.ncrit)

        # ========== NACA POINTS DEDICATED TEST ==========
        naca_grp = f.create_group('naca_points')
        # Test various NACA airfoils
        naca_codes = ['0012', '2412', '4415', '0006']
        for code in naca_codes:
            M_naca = mfoil.mfoil(naca=code, npanel=100)
            grp = naca_grp.create_group(code)
            grp.create_dataset('xpoint', data=M_naca.geom.xpoint)
            grp.create_dataset('npoint', data=M_naca.geom.npoint)

        # ========== SET COORDS DEDICATED TEST ==========
        setcoords_grp = f.create_group('set_coords')
        # Test SetCoords by creating a foil and checking the result
        M_coords = mfoil.mfoil(naca='0012', npanel=130)
        setcoords_grp.create_dataset('x', data=M_coords.foil.x)
        setcoords_grp.create_dataset('s', data=M_coords.foil.s)
        setcoords_grp.create_dataset('t', data=M_coords.foil.t)
        setcoords_grp.create_dataset('N', data=M_coords.foil.N)
        setcoords_grp.create_dataset('xpoint', data=M_coords.geom.xpoint)

        # ========== SPLINE 1D CREATE DEDICATED TEST ==========
        spline1d_grp = f.create_group('spline1d')
        # Create 1D spline from arc length vs x coordinate
        from scipy.interpolate import CubicSpline
        s_vals = M_coords.foil.s
        x_vals = M_coords.foil.x[0, :]  # x coordinates
        y_vals = M_coords.foil.x[1, :]  # y coordinates

        # Create spline
        spline_x = CubicSpline(s_vals, x_vals)
        spline_y = CubicSpline(s_vals, y_vals)

        # Save spline data
        spline1d_grp.create_dataset('s', data=s_vals)
        spline1d_grp.create_dataset('x_vals', data=x_vals)
        spline1d_grp.create_dataset('y_vals', data=y_vals)
        spline1d_grp.create_dataset('x_coeffs', data=spline_x.c)
        spline1d_grp.create_dataset('y_coeffs', data=spline_y.c)

        # Evaluate at test points
        s_test = np.linspace(0, s_vals[-1], 50)
        x_eval = spline_x(s_test)
        y_eval = spline_y(s_test)
        dx_eval = spline_x(s_test, 1)  # First derivative
        dy_eval = spline_y(s_test, 1)
        d2x_eval = spline_x(s_test, 2)  # Second derivative
        d2y_eval = spline_y(s_test, 2)

        spline1d_grp.create_dataset('s_test', data=s_test)
        spline1d_grp.create_dataset('x_eval', data=x_eval)
        spline1d_grp.create_dataset('y_eval', data=y_eval)
        spline1d_grp.create_dataset('dx_eval', data=dx_eval)
        spline1d_grp.create_dataset('dy_eval', data=dy_eval)
        spline1d_grp.create_dataset('d2x_eval', data=d2x_eval)
        spline1d_grp.create_dataset('d2y_eval', data=d2y_eval)

        # ========== INIT THERMO COMPLETE REFERENCE ==========
        # Save all thermo parameters for comprehensive InitThermo test
        init_thermo = f.create_group('init_thermo')

        # Test case 1: Standard conditions (from M above which was set up with alpha=0, Re=1e6)
        init_thermo.create_dataset('Vinf', data=M.param.Vinf)
        init_thermo.create_dataset('Minf', data=M.param.Minf)
        init_thermo.create_dataset('muinf', data=M.param.muinf)
        init_thermo.create_dataset('mu0', data=M.param.mu0)
        init_thermo.create_dataset('rho0', data=M.param.rho0)
        init_thermo.create_dataset('H0', data=M.param.H0)
        init_thermo.create_dataset('KTb', data=M.param.KTb)
        init_thermo.create_dataset('KTl', data=M.param.KTl)
        init_thermo.create_dataset('cps', data=M.param.cps)
        init_thermo.create_dataset('gam', data=M.param.gam)
        init_thermo.create_dataset('Re', data=M.oper.Re)
        init_thermo.create_dataset('Ma', data=M.oper.Ma)

        # Test case 2: Compressible flow
        M_comp = mfoil.mfoil(naca='0012', npanel=130)
        M_comp.setoper(alpha=0.0, Re=1e6, Ma=0.5)
        M_comp.param.verb = 0
        mfoil.init_thermo(M_comp)

        comp_grp = init_thermo.create_group('compressible')
        comp_grp.create_dataset('Ma', data=M_comp.oper.Ma)
        comp_grp.create_dataset('KTb', data=M_comp.param.KTb)
        comp_grp.create_dataset('KTl', data=M_comp.param.KTl)
        comp_grp.create_dataset('cps', data=M_comp.param.cps)

    print(f"Reference data saved to {filename}")
    return M


def main():
    print("Extracting mfoil reference data to HDF5...")
    output_file = '/Users/clupp/tools/analysis/ProjectTrysil/cppfoil/tests/mfoil_reference.h5'
    M = extract_all_reference_data(output_file)

    # Print summary
    with h5py.File(output_file, 'r') as f:
        print(f"\nSummary:")
        print(f"  Geometry: N={f['geometry/N'][()]}")
        print(f"  Panel influence (source) tests: {f['panel_influence'].attrs['num_tests']}")
        print(f"  Panel influence (vortex) tests: {f['panel_vortex'].attrs['num_tests']}")
        print(f"  Panel info tests: {f['panel_info'].attrs['num_tests']}")
        print(f"  Stagnation: Istag={f['stagnation/Istag'][:]}")
        print(f"  Wake: N={f['wake/N'][()]}")
        print(f"  Surfaces: lower={len(f['surfaces/Is_lower'][:])}, upper={len(f['surfaces/Is_upper'][:])}, wake={len(f['surfaces/Is_wake'][:])}")
        print(f"  BL property tests: {len(list(f['bl_properties'].keys()))}")
        print(f"  Residual tests: {f['residual_station'].attrs['num_tests']}")
        print(f"  Initial |R|: {f['glob_sys/R_norm'][()]:.6f}")
        print(f"  Inviscid velocity tests: {f['inviscid_velocity'].attrs['num_tests']}")
        print(f"  Spline: quadseg_w sum={f['spline/quadseg_w'][:].sum():.6f}")
        print(f"  Forces: cl={f['forces/cl'][()]:.6f}, cd={f['forces/cd'][()]:.6f}")


if __name__ == '__main__':
    main()
