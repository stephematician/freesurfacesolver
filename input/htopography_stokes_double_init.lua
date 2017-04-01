-- dofile('input/htopography_find_stokes_double.lua')

io.write('Computation for the following desired parameters:\n')
io.write('       B : ', B, '\n')
io.write('       L : ', L, '\n')
if not TYPE4_FROUDE then
  io.write('  FROUDE : ', FROUDE, '\n\n')
else
  io.write('       A : ', TYPE4_A, '\n\n')
end

io.write('1. Find a low-resolution stokes double wave solution.\n')
io.flush()

if TYPE4_A and TYPE4_FROUDE then
lo_res_dbl,
free_lo_res_dbl,
continued_lo_res_dbl = find_stokes_double(TYPE4_A,
                                                B,
                                                L,
                                     STOKES_U_MAX,
                                          PHI_LIM,
                                         N_LO_RES,
                                 CLUSTER_DPHI_VAL,
                                     TYPE4_FROUDE)

else
lo_res_dbl,
free_lo_res_dbl,
continued_lo_res_dbl = find_stokes_double(DOUBLE_A,
                                                 B,
                                                 L,
                                      STOKES_U_MIN,
                                           PHI_LIM,
                                          N_LO_RES,
                                  CLUSTER_DPHI_VAL)

end
io.write('----------------\n')
io.write('2. Find a high-resolution stokes double wave solution.\n')
io.flush()

     lo_res_beta_s,
lo_res_theta_s_mid,
  lo_res_tau_s_mid,
    lo_res_x_s_mid,
      lo_res_x_lhs,
  lo_res_eta_s_mid,
    lo_res_eta_lhs = htopography_residuals.compute_bttxe_hs(lo_res_dbl)

lo_res_phi_s = integrals_sbox.phi(lo_res_beta_s,
                            lo_res_dbl.BETA_SUB,
                             lo_res_dbl.PHI_SUB,
                            lo_res_dbl.DPHI_SUB,
                          lo_res_dbl.HOMOTOPY_S)

io.write('  Interpolating onto high resolution grid ... ')
io.flush()

hi_res_dbl           = htopography_grid.deepcopy(lo_res_dbl)
free_hi_res_dbl      = htopography_grid.deepcopy(free_lo_res_dbl)
continued_hi_res_dbl = htopography_grid.deepcopy(continued_lo_res_dbl)

 hi_res_dbl.THETA_S,
hi_res_dbl.BETA_SUB,
   hi_res_dbl.N_SUB,
 hi_res_dbl.PHI_SUB,
hi_res_dbl.DPHI_SUB = grid_functions.convert_grid_hs(lo_res_phi_s,
                                               lo_res_dbl.THETA_S,
                                                 lo_res_tau_s_mid,
                                                         N_HI_RES,
                                                 CLUSTER_DPHI_VAL,
                                            lo_res_dbl.HOMOTOPY_S,
                                               integrals_sbox.phi)

io.write('done.\n')
io.flush()

free_hi_res_dbl.THETA_S = {}
for i = 1, N_HI_RES do
  free_hi_res_dbl.THETA_S[i] = true
end
free_hi_res_dbl.THETA_S[1] = false

continued_hi_res_dbl.THETA_S = {}
for i = 1, N_HI_RES do
  continued_hi_res_dbl.THETA_S[i] = false
end

if TYPE4_A and TYPE4_FROUDE then
  free_hi_res_dbl.A = false
  free_hi_res_dbl.U_MIN = true

  continued_hi_res_dbl.A = false
  continued_hi_res_dbl.U_MIN = true
end

io.write('  Creating new (high resolution) surface object ... ')
io.flush()

freesurface = surface.new(     hi_res_dbl,
                          free_hi_res_dbl,
                     continued_hi_res_dbl,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)

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
   -- htopography_residuals.residual_hs(freesurface.unknowns, {freesurface})

   error('Could not find new higher-resolution surface.')
end

io.write('  Parameters for high resolution stokes double wave,\n')
io.write('       N :')
for k,v in pairs(hi_res_dbl.N_SUB) do
   io.write(' ' .. tostring(v))
end
io.write(',\n')
io.write('       A : ', hi_res_dbl.A, ',\n')
io.write('       B : ', hi_res_dbl.B, ',\n')
io.write('       D : ', hi_res_dbl.D, ',\n')
io.write('  FROUDE : ', hi_res_dbl.FROUDE, ',\n')
io.write('  LAMBDA : ', hi_res_dbl.LAMBDA, ',\n')
io.write('    ETA0 : ', hi_res_dbl.ETA0, ',\n')
io.write('   PHI_C : ', hi_res_dbl.PHI_C, '.\n\n')
