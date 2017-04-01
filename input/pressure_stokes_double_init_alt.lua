dofile('input/pressure_find_doublelim_from_unif.lua')

io.write('Computation for the following desired parameters\n')
io.write('      U(0) : ', STOKES_U_MIN, '\n')
io.write('         B : ', B, '\n')

io.write('1. Find a low-resolution double hump solution with small u.\n')
io.flush()

if TYPE4_A and TYPE4_FROUDE then
lo_res_dbl,
free_lo_res_dbl,
continued_lo_res_dbl = find_pert_to_double(B,
                                STOKES_U_MAX,
                                     PHI_LIM,
                                    N_LO_RES,
                            CLUSTER_DPHI_VAL,
                                     TYPE4_A,
                                TYPE4_FROUDE)

else
lo_res_dbl,
free_lo_res_dbl,
continued_lo_res_dbl = find_pert_to_double(B,
                                STOKES_U_MIN,
                                     PHI_LIM,
                                    N_LO_RES,
                            CLUSTER_DPHI_VAL)

end

io.write('----------------\n')
io.write('2. Find a high-resolution unforced solitary wave solution.\n')
io.flush()

-- find_pert_to_sol has incorrect free/continuation variables - probably need to fix this up later!

io.write('----------------\n')
io.write('Refining surface\n')
io.flush()

   lo_res_beta_s,
  lo_res_u_s_mid,
  lo_res_v_s_mid,
  lo_res_x_s_mid,
    lo_res_x_lhs,
lo_res_eta_s_mid,
  lo_res_eta_lhs = pressure_residuals.compute_buvxe_hs(lo_res_dbl)

lo_res_phi_s = integrals_pressure.phi(lo_res_beta_s,
                                lo_res_dbl.BETA_SUB,
                                 lo_res_dbl.PHI_SUB,
                                lo_res_dbl.DPHI_SUB,
                              lo_res_dbl.HOMOTOPY_S)

io.write('  Interpolating onto high resolution grid ... ')
io.flush()

-- pressure grid needs deepcopy function?
hi_res_dbl           = pressure_grid.deepcopy(lo_res_dbl)
free_hi_res_dbl      = pressure_grid.deepcopy(free_lo_res_dbl)
continued_hi_res_dbl = pressure_grid.deepcopy(continued_lo_res_dbl)

     hi_res_dbl.V_S,
hi_res_dbl.BETA_SUB,
   hi_res_dbl.N_SUB,
 hi_res_dbl.PHI_SUB,
hi_res_dbl.DPHI_SUB = grid_functions.convert_grid_uv_hs(lo_res_phi_s,
                                                      lo_res_dbl.V_S,
                                                      lo_res_u_s_mid,
                                                            N_HI_RES,
                                                    CLUSTER_DPHI_VAL,
                                               lo_res_dbl.HOMOTOPY_S,
                                              integrals_pressure.phi)
io.write('done.\n')
io.flush()

free_hi_res_dbl.V_S = {}
for i = 1, N_HI_RES do
  free_hi_res_dbl.V_S[i] = true
end
free_hi_res_dbl.V_S[1] = false

continued_hi_res_dbl.THETA_S = {}
for i = 1, N_HI_RES do
  continued_hi_res_dbl.THETA_S[i] = false
end

free_hi_res_dbl.A = false
free_hi_res_dbl.U_MIN = true

continued_hi_res_dbl.A = false
continued_hi_res_dbl.U_MIN = true

io.write('  Creating new (high resolution) surface object ... ')
io.flush()

freesurface = surface.new(    hi_res_dbl,
                         free_hi_res_dbl,
                    continued_hi_res_dbl,
          pressure_residuals.residual_hs,
    pressure_residuals.compute_output_hs)

io.write('done.\n')

DIRECTION = continuator.BACKWARD

io.write('  Use Newton\'s method to solve on the new grid ... ')
io.flush()
res = freesurface:progress(0, DIRECTION, freesurface)
io.write('done.\n')

hi_res_dbl = surface.unpack_vec(freesurface.unknowns,
                                 freesurface.initial,
                                 freesurface.is_free,
                            freesurface.is_continued,
                               freesurface.ind_table)

if res ~= continuator.CONVERGED then
   -- pressure_residuals.residual_hs(freesurface.unknowns, {freesurface})

   error('Could not find new higher-resolution surface.')
end

io.write('  Parameters for initial higher resolution surface,\n')
io.write('             n : ')
for k,v in pairs(hi_res_dbl.N_SUB) do
   io.write(' ' .. tostring(v))
end
io.write(',\n')
io.write('       A : ', hi_res_dbl.A, ',\n')
io.write('       B : ', hi_res_dbl.B, ',\n')
io.write('       D : ', hi_res_dbl.D, ',\n')
io.write('  LAMBDA : ', hi_res_dbl.LAMBDA, ',\n')
io.write('    ETA0 : ', hi_res_dbl.ETA0, '.\n')
io.write('       F : ', hi_res_dbl.FROUDE, '.\n')
io.write('   U_MIN : ', hi_res_dbl.U_MIN, '.\n\n')
