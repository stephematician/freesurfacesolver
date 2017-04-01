require "surface"
require "pressure_residuals"
require "pressure_grid"

dofile('input/pressure_find_pert_to_unif.lua')

N = 800
N_ROUGH = 400

desired_A = -0.20

B                = 2.80
FROUDE           = 1.262
PHI_LIM          = 11
CLUSTER_DPHI_VAL = 0.0025

MAX_AREA = 2*B

STOKES_U_MIN = 0.05

display_initial = false

io.write('Computation for the following desired parameters\n')
io.write('        a : ', desired_A, '\n')
io.write('        b : ', B, '\n')
io.write('   froude : ', FROUDE, '\n\n')

io.write('Creating new surface.\n')
io.flush()

initial_surf,
is_free_surf,
is_continued_surf = find_pert_to_unif(FROUDE,
                                     PHI_LIM,
                                           B,
                                     N_ROUGH,
                            CLUSTER_DPHI_VAL,
                                   desired_A)

io.write('----------------\n')
io.write('Refining surface\n')
io.flush()

ibeta_s, iu_s_mid, iv_s_mid,
         ix_s_mid,   ix_lhs,
       ieta_s_mid, ieta_lhs = pressure_residuals.compute_buvxe_hs(initial_surf)

iphi_s  = integrals_pressure.phi(ibeta_s,
                   initial_surf.BETA_SUB,
                    initial_surf.PHI_SUB,
                   initial_surf.DPHI_SUB,
                 initial_surf.HOMOTOPY_S)

     initial_surf.V_S,
initial_surf.BETA_SUB,
   initial_surf.N_SUB,
 initial_surf.PHI_SUB,
initial_surf.DPHI_SUB = grid_functions.convert_grid_uv_hs(iphi_s,
                                                initial_surf.V_S,
                                                        iu_s_mid,
                                                               N,
                                                CLUSTER_DPHI_VAL,
                                         initial_surf.HOMOTOPY_S,
                                          integrals_pressure.phi)

is_free_surf.V_S = {}
for i = 1, N do
  is_free_surf.V_S[i] = true
end
is_free_surf.V_S[1] = false

is_continued_surf.V_S = {}
for i = 1, N do
  is_continued_surf.V_S[i] = false
end

io.write('  Creating new surface object ... ')
io.flush()

freesurface = surface.new(   initial_surf,
                             is_free_surf,
                        is_continued_surf,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

io.write('done.\n')

DIRECTION = continuator.BACKWARD

io.write('  Attempting to converge upon new higher-resolution surface...')
io.flush()
res = freesurface:progress(0, DIRECTION, freesurface)
io.write('done.\n')

initial_surf = surface.unpack_vec(freesurface.unknowns,
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
for k,v in pairs(initial_surf.N_SUB) do
   io.write(tostring(v) .. ' ')
end
io.write(',\n')
io.write('       A : ', initial_surf.A, ',\n')
io.write('       B : ', initial_surf.B, ',\n')
io.write('       D : ', initial_surf.D, ',\n')
io.write('  lambda : ', initial_surf.LAMBDA, ',\n')
io.write('    eta0 : ', initial_surf.ETA0, '.\n\n')


