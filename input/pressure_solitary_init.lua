-- dofile('input/pressure_find_pert_to_unif.lua')
dofile('input/pressure_find_sollim_from_unif.lua')

io.write('Computation for the following desired parameters\n')
io.write('     u(0) : ', STOKES_U_MIN, '\n')
io.write('        b : ', B, '\n')
io.write('   froude : ', FROUDE, '\n\n')


io.write('Creating new surface.\n')
io.flush()

initial_sol,
is_free_sol,
is_continued_sol = find_pert_to_sol(PHI_LIM,
                                          B,
                                    N_ROUGH,
                           CLUSTER_DPHI_VAL,
                               STOKES_U_MIN,
                                     FROUDE)

-- find_pert_to_sol has incorrect free variables
is_free_sol.U_MIN = false
is_free_sol.FROUDE = true

io.write('----------------\n')
io.write('Refining surface\n')
io.flush()

ibeta_s, iu_s_mid, iv_s_mid,
         ix_s_mid,   ix_lhs,
       ieta_s_mid, ieta_lhs = pressure_residuals.compute_buvxe_hs(initial_sol)

iphi_s  = integrals_pressure.phi(ibeta_s,
                    initial_sol.BETA_SUB,
                     initial_sol.PHI_SUB,
                    initial_sol.DPHI_SUB,
                  initial_sol.HOMOTOPY_S)

     initial_sol.V_S,
initial_sol.BETA_SUB,
   initial_sol.N_SUB,
 initial_sol.PHI_SUB,
initial_sol.DPHI_SUB = grid_functions.convert_grid_uv_hs(iphi_s,
                                                initial_sol.V_S,
                                                       iu_s_mid,
                                                              N,
                                               CLUSTER_DPHI_VAL,
                                         initial_sol.HOMOTOPY_S,
                                         integrals_pressure.phi)

is_free_sol.V_S = {}
for i = 1, N do
  is_free_sol.V_S[i] = true
end
is_free_sol.V_S[1] = false

is_continued_sol.V_S = {}
for i = 1, N do
  is_continued_sol.V_S[i] = false
end

io.write('  Creating new surface object ... ')
io.flush()

freesurface = surface.new(    initial_sol,
                              is_free_sol,
                         is_continued_sol,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

io.write('done.\n')

DIRECTION = continuator.BACKWARD

io.write('  Attempting to converge upon new higher-resolution surface...')
io.flush()
res = freesurface:progress(0, DIRECTION, freesurface)
io.write('done.\n')

initial_sol = surface.unpack_vec(freesurface.unknowns,
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
for k,v in pairs(initial_sol.N_SUB) do
   io.write(tostring(v) .. ' ')
end
io.write(',\n')
io.write('       A : ', initial_sol.A, ',\n')
io.write('       B : ', initial_sol.B, ',\n')
io.write('       D : ', initial_sol.D, ',\n')
io.write('  lambda : ', initial_sol.LAMBDA, ',\n')
io.write('    eta0 : ', initial_sol.ETA0, '.\n\n')

