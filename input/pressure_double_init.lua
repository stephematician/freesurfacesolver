dofile('input/pressure_find_sollim_from_unif.lua')

io.write('Computation for the following desired parameters:\n')
io.write('       B : ', B, '\n')
io.write('  FROUDE : ', FROUDE, '\n\n')

io.write('1. Find a low-resolution perturbation to solitary wave solution.\n')
io.flush()

lo_res_sol,
free_lo_res_sol,
continued_lo_res_sol = find_pert_to_sol(B,
                                   FROUDE,
                             STOKES_U_MIN,
                                  PHI_LIM,
                                 N_LO_RES,
		         CLUSTER_DPHI_VAL)


io.write('----------------\n')
io.write('2. Find a high-resolution  perturbation to double hump solution.\n')
io.flush()

   lo_res_beta_s,
  lo_res_u_s_mid,
  lo_res_v_s_mid,
  lo_res_x_s_mid,
    lo_res_x_lhs,
lo_res_eta_s_mid,
  lo_res_eta_lhs = pressure_residuals.compute_buvxe_hs(lo_res_sol)

lo_res_phi_s  = integrals_pressure.phi(lo_res_beta_s,
                                 lo_res_sol.BETA_SUB,
                                  lo_res_sol.PHI_SUB,
                                 lo_res_sol.DPHI_SUB,
                               lo_res_sol.HOMOTOPY_S)

io.write('  Interpolating onto high resolution grid ... ')
io.flush()

hi_res_sol           = pressure_grid.deepcopy(lo_res_sol)
free_hi_res_sol      = pressure_grid.deepcopy(free_lo_res_sol)
continued_hi_res_sol = pressure_grid.deepcopy(continued_lo_res_sol)

     hi_res_sol.V_S,
hi_res_sol.BETA_SUB,
   hi_res_sol.N_SUB,
 hi_res_sol.PHI_SUB,
hi_res_sol.DPHI_SUB = grid_functions.convert_grid_uv_hs(lo_res_phi_s,
                                                      lo_res_sol.V_S,
                                                      lo_res_u_s_mid,
                                                            N_HI_RES,
                                                    CLUSTER_DPHI_VAL,
                                               lo_res_sol.HOMOTOPY_S,
                                              integrals_pressure.phi)

io.write('done.\n')
io.flush()

free_hi_res_sol.V_S = {}
for i = 1, N_HI_RES do
  free_hi_res_sol.V_S[i] = true
end
free_hi_res_sol.V_S[1] = false

continued_hi_res_sol.V_S = {}
for i = 1, N_HI_RES do
  continued_hi_res_sol.V_S[i] = false
end


io.write('  Creating new (high resolution) surface object ... ')
io.flush()

freesurface = surface.new(  hi_res_sol,
                       free_hi_res_sol,
                  continued_hi_res_sol,
        pressure_residuals.residual_hs,
  pressure_residuals.compute_output_hs)

io.write('done.\n')

DIRECTION = continuator.BACKWARD

io.write('  Use Newton\'s method to solve on the new grid ... ')
io.flush()
res = freesurface:progress(0, DIRECTION, freesurface)
io.write('done.\n')

hi_res_sol = surface.unpack_vec(freesurface.unknowns,
                                 freesurface.initial,
                                 freesurface.is_free,
                            freesurface.is_continued,
                               freesurface.ind_table)

if res ~= continuator.CONVERGED then
   -- pressure_residuals.residual_hs(freesurface.unknowns, {freesurface})

   error('Could not find new higher-resolution surface.')
end


io.write('  Parameters for high resolution perturbation to double hump,\n')
io.write('       N :')
for k,v in pairs(hi_res_sol.N_SUB) do
   io.write(' ' .. tostring(v))
end
io.write(',\n')
io.write('       A : ', hi_res_sol.A, ',\n')
io.write('       B : ', hi_res_sol.B, ',\n')
io.write('       D : ', hi_res_sol.D, ',\n')
io.write('  LAMBDA : ', hi_res_sol.LAMBDA, ',\n')
io.write('    ETA0 : ', hi_res_sol.ETA0, ',\n')
io.write('       F : ', hi_res_sol.FROUDE, ',\n')
io.write('   U_MIN : ', hi_res_sol.U_MIN, '.\n\n')

