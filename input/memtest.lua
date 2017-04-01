--[[
Run a basic system and solve for it, then run one continuation step to test memory leaks

Box is defined via the use of an approximation to a step function given by
the hyperbolic tangent.

homotopy available, symmetric formulation, with dbeta dependent
upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 21/06/2014
--]]

require "surface"
require "htopography_residuals"
require "htopography_grid"

N_LO_RES = 40

L       = 2.00
B       = 10.00
FROUDE  = 1.10
PHI_LIM = 10

CLUSTER_DPHI_VAL = 0.01


  local approx_LAMBDA = math.sqrt((-5/4) + 
                                  math.sqrt((25/16) + 
                                            (15/2)*(FROUDE*FROUDE - 1)))
  io.write('\nSolved approximate equation; LAMBDA = ', approx_LAMBDA, ' + h.o.t. \n\n')
  io.flush()

  local DPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- A and D are 'magic' numbers that get us started near a uniform stream
  local initial_fs = { A = 1e-5,
                       B = B,
                       L = L,
                       FROUDE = FROUDE,
                       LAMBDA = approx_LAMBDA,
                       PHI_C = L/2.0,
                       U_MIN = 1.0,
                       D = 0.0005,
                       THETA_S = {},
                       BETA_SUB = { 1, N_LO_RES },
                       PHI_SUB = {0, PHI_LIM},
                       DPHI_SUB = {DPHI_LIM/N_LO_RES, DPHI_LIM/N_LO_RES},
                       N_SUB = { N_LO_RES-1 },
                       CLUSTER_VAL = CLUSTER_DPHI_VAL,
                       HOMOTOPY_S = 1.0,
                       GAMMA = 25.0, -- see note??
                       ETA0 = 0
                     }

  for j = 1, N_LO_RES do
    initial_fs.THETA_S[j] = 0
  end

  local free_lo_res_unif = { A = true,
                             B = false,
                             L = false,
                             FROUDE = false,
                             LAMBDA = true,
                             PHI_C = true,
                             U_MIN = true,
                             D = true,
                             THETA_S = {},
                             BETA_SUB = {false, false},
                             PHI_SUB = {false, false},
                             DPHI_SUB = {true, true},
                             N_SUB = {false},
                             CLUSTER_VAL = false,
                             HOMOTOPY_S = false,
                             GAMMA = false,
                             ETA0 = true
                           }

  for i = 1, N_LO_RES do
    free_lo_res_unif.THETA_S[i] = true
  end
  free_lo_res_unif.THETA_S[1] = false

  --[[ The next table 'is_continued' specifies which variables are to be continued
     from. At the moment ONE and ONLY ONE value may be set to true, the rest must
     be false --]]

  local continued_lo_res_unif = { A = true,
                                  B = false,
                                  L = false,
                                  FROUDE = false,
                                  LAMBDA = false,
                                  PHI_C = false,
                                  U_MIN = false,
                                  D = false,
                                  THETA_S = {},
                                  BETA_SUB = {false, false},
                                  PHI_SUB = {false, false},
                                  DPHI_SUB = {false, false},
                                  N_SUB = {false},
                                  CLUSTER_VAL = false,
                                  HOMOTOPY_S = false,
                                  GAMMA = false,
                                  ETA0 = false
                                }
  for i = 1, N_LO_RES do
   continued_lo_res_unif.THETA_S[i] = false
  end

  io.write('  Creating new surface object ... ')
  io.flush()
  local freesurface = surface.new(initial_fs,
                            free_lo_res_unif,
                       continued_lo_res_unif,
           htopography_residuals.residual_hs,
     htopography_residuals.compute_output_hs)
  io.write('done.\n')

  local DIRECTION = continuator.BACKWARD

  io.write('  Attempting to find initial near uniform stream solution ... ')
  io.flush()
  local res = freesurface:progress(0, DIRECTION, freesurface)
  io.write('done.\n')
  io.flush()
  io.write(res)
  local lo_res_unif = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)

-- Done with new surface?

freesurface = surface.new(    lo_res_unif,
                         free_lo_res_unif,
                    continued_lo_res_unif,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)

freesurface:set_length_vars({'A', 'ETA0'})

ds = 0.01

stepsize = 0.01
max_stepsize = 2.5

freesurface:set_max_stepsize(max_stepsize)
freesurface:set_stepsize(stepsize)

local progress_start = os.clock()
--  res, msg = freesurface:progress_interp(ds, DIRECTION, freesurface)
local progress_finish = os.clock()

--io.write(msg)
io.write('    Continuation time = ' .. tostring(progress_finish - progress_start) ..'\n')
io.flush()

local regrid_start = os.clock()
--freesurface, direction = htopography_grid.regrid_hs(freesurface)
local regrid_finish = os.clock()
io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
io.flush()

collectgarbage('collect')

io.write('Finished continuation and regrid step.\n')
io.flush()
