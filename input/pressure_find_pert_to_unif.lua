find_pert_to_unif = function(finalA, inB, inFROUDE, inPHI_LIM, inN, inCLUSTER_VAL)
  --[[ \todo Write some documentation here ]]--



  io.write('\n--------------------------------\n')
  io.write('Running script FIND_PERT_TO_UNIF\n')
  io.write('--------------------------------\n')
  io.write('Begin with a near uniform-stream solution with gaussian pressure disturbance.\n')
  io.write('Pressure parameters, B = ', inB, ', and ')
  io.write('continue to nearby solutions with A = ', finalA, '.\n')

  local approx_LAMBDA = math.sqrt((-5/4) + 
                                  math.sqrt((25/16) + 
                                            (15/2)*(inFROUDE*inFROUDE - 1)))
  io.write('\nSolved approximate equation; LAMBDA = ', approx_LAMBDA, ' + h.o.t. \n\n')
  io.flush()

  local inDPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- A and D are 'magic' numbers that get us started near a uniform stream
  local initial_fs = { A = 1e-5,
                       B = inB,
                       FROUDE = inFROUDE,
                       LAMBDA = approx_LAMBDA,
                       U_MIN = 1.0,
                       D = 0.0005,
                       V_S = {},
                       BETA_SUB = { 1, inN },
                       PHI_SUB = {0, inPHI_LIM},
                       DPHI_SUB = {inDPHI_LIM/inN, inDPHI_LIM/inN},
                       N_SUB = { inN-1 },
                       CLUSTER_VAL = inCLUSTER_VAL,
                       HOMOTOPY_S = 1.0,
                       ETA0 = 0
                     }

  for j = 1, inN do
    initial_fs.V_S[j] = 0
  end

  local is_free_fs = { A = true,
                       B = false,
                       FROUDE = false,
                       LAMBDA = true,
                       U_MIN = true,
                       D = true,
                       V_S = {},
                       BETA_SUB = {false, false},
                       PHI_SUB = {false, false},
                       DPHI_SUB = {true, true},
                       N_SUB = {false},
                       CLUSTER_VAL = false,
                       HOMOTOPY_S = false,
                       ETA0 = true
                     }

  for i = 1, inN do
    is_free_fs.V_S[i] = true
  end
  is_free_fs.V_S[1] = false

  --[[ The next table 'is_continued' specifies which variables are to be continued
     from. At the moment ONE and ONLY ONE value may be set to true, the rest must
     be false ]]--

  local is_continued_fs = { A = true,
                            B = false,
                            FROUDE = false,
                            LAMBDA = false,
                            U_MIN = false,
                            D = false,
                            V_S = {},
                            BETA_SUB = {false, false},
                            PHI_SUB = {false, false},
                            DPHI_SUB = {false, false},
                            N_SUB = {false},
                            CLUSTER_VAL = false,
                            HOMOTOPY_S = false,
                            ETA0 = false
                          }
  for i = 1, inN do
    is_continued_fs.V_S[i] = false
  end

  io.write('  Creating new surface object ... ')
  io.flush()

  local freesurface = surface.new(initial_fs,
                                  is_free_fs,
                             is_continued_fs,
              pressure_residuals.residual_hs,
        pressure_residuals.compute_output_hs)

  io.write('done.\n')
  io.flush()

  local DIRECTION = continuator.BACKWARD

  io.write('  Attempting to find initial near uniform stream solution ... ')
  io.flush()
  local res = freesurface:progress(0, DIRECTION, freesurface)
  io.write('done.\n')
  io.flush()

  local initial_unif = surface.unpack_vec(freesurface.unknowns,
                                           freesurface.initial,
                                           freesurface.is_free,
                                      freesurface.is_continued,
                                         freesurface.ind_table)

  io.write('  Parameters for initial near uniform stream solution,\n')
  io.write('       n : ')
  for k,v in pairs(initial_unif.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('       A : ', initial_unif.A, ',\n')
  io.write('       B : ', initial_unif.B, ',\n')
  io.write('       D : ', initial_unif.D, ',\n')
  io.write('  lambda : ', initial_unif.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_unif.ETA0, '.\n\n')

  if res ~= continuator.CONVERGED then
     print_out = true
     free_surface.residual(freesurface.unknowns, {freesurface})

     error('Could not find initial guassian topography solution.')
  end

  io.write('Begin continuation with respect to A.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 100
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.6
  local started_converge = false
  local ds = 0.03
  local max_stepsize = 0.3
  local initial_step = true

  local res = nil
  local msg = nil

  freesurface:set_max_stepsize(max_stepsize)
  freesurface:set_stepsize(stepsize)

  for i = 1, num_continuate do
    local outstr = 'Iteration : ' .. tostring(i-1)
    local filstr = '-'
    io.write(outstr .. '\n')
    io.write(filstr:rep(outstr:len()) .. '\n')

    local csurf = surface.unpack_vec(freesurface.unknowns,
                                      freesurface.initial,
                                      freesurface.is_free,
                                 freesurface.is_continued,
                                    freesurface.ind_table)

    io.write('  Parameter summary :\n')  
    io.write('       A : ', csurf.A, '\n')
    io.write('    ETA0 : ', csurf.ETA0, '\n')

    if (not initial_step) and
       (math.abs(csurf.A) > (math.abs(finalA) - ds)) then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
      end
      res, msg =  freesurface:progress_approach(finalA, rough_converge_rate, freesurface)
      io.write(msg)

    else

      res, msg = freesurface:progress_interp(ds, direction, freesurface)
      io.write(msg)

    end

    if res == continuator.BIFURCATION then
      io.write('    Bifurcation detected, attempting to switch.\n')

      if direction == continuator.FORWARD then
        direction = continuator.BACKWARD
      else
        direction = continuator.FORWARD
      end

    end

    if res == continuator.ERROR then
      error('Unexpected error detected in continuator:progress().')
      break
    end

    if res == continuator.DIVERGED then
      error('Unexpected divergence detected in continuator:progress().')
      break
    end

    if initial_step then
      local nsurf = surface.unpack_vec(freesurface.unknowns,
                                        freesurface.initial,
                                        freesurface.is_free,
                                    freesurface.is_continued,
                                       freesurface.ind_table)
      local lsgn = 1
      if finalA < 0 then
        lsgn = -1
      end
      if lsgn*(csurf.A) > lsgn*(nsurf.A) then
        freesurface:reset_interp()
        io.write('    Started in wrong direction, attempting to switch.\n')
        if direction == continuator.FORWARD then
          direction = continuator.BACKWARD
        else
          direction = continuator.FORWARD
        end
      else
        initial_step = false
      end 
    end

    if started_converge and (i - started_converge_i) > 7 then
      io.write('  Finished.\n')
      break
    end

  end

  local initial_pert = surface.unpack_vec(freesurface.unknowns,
                                           freesurface.initial,
                                           freesurface.is_free,
                                      freesurface.is_continued,
                                         freesurface.ind_table)

  initial_pert.A = finalA

  return initial_pert, is_free_fs, is_continued_fs

end