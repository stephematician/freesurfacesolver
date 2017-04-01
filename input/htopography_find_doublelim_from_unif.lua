pert_to_double_tables = function(inN, inN_SUB_length)

  local initial_fs = { A = 0,
                       B = 0,
                       L = 0,
                       FROUDE = 0,
                       LAMBDA = 0,
                       PHI_C = 0,
                       U_MIN = 0,
                       D = 0,
                       THETA_S = {},
                       BETA_SUB = {},
                       PHI_SUB = {},
                       DPHI_SUB = {},
                       N_SUB = {},
                       CLUSTER_VAL = 0,
                       HOMOTOPY_S = 0,
                       GAMMA = 0,
                       ETA0 = 0
                     }

  for j = 1, inN do
    initial_fs.THETA_S[j] = nil
  end

  local is_free_fs = { A = true,
                       B = false,
                       L = false,
                       FROUDE = false,
                       LAMBDA = true,
                       PHI_C = true,
                       U_MIN = true,
                       D = true,
                       THETA_S = {},
                       BETA_SUB = {},
                       PHI_SUB = {},
                       DPHI_SUB = {},
                       N_SUB = {},
                       CLUSTER_VAL = false,
                       HOMOTOPY_S = false,
                       GAMMA = false,
                       ETA0 = true
                     }

  for i = 1, inN do
    is_free_fs.THETA_S[i] = true
  end
  is_free_fs.THETA_S[1] = false

  --[[ The next table 'is_continued' specifies which variables are to be continued
     from. At the moment ONE and ONLY ONE value may be set to true, the rest must
     be false --]]

  local is_continued_fs = { A = true,
                            B = false,
                            L = false,
                            FROUDE = false,
                            LAMBDA = false,
                            PHI_C = false,
                            U_MIN = false,
                            D = false,
                            THETA_S = {},
                            BETA_SUB = {},
                            PHI_SUB = {},
                            DPHI_SUB = {},
                            N_SUB = {},
                            CLUSTER_VAL = false,
                            HOMOTOPY_S = false,
                            GAMMA = false,
                            ETA0 = false
                          }
  for i = 1, inN do
    is_continued_fs.THETA_S[i] = false
  end

  is_free_fs.PHI_SUB = { false }
  is_free_fs.BETA_SUB = { false }
  is_free_fs.DPHI_SUB = { true }
  is_free_fs.N_SUB = { false }
  is_continued_fs.PHI_SUB = { false }
  is_continued_fs.BETA_SUB = { false }
  is_continued_fs.DPHI_SUB = { false }
  is_continued_fs.N_SUB = { false }

  for i = 1, inN_SUB_length-1 do
    is_free_fs.BETA_SUB[i+1] = false
    is_free_fs.N_SUB[i+1] = false
    is_free_fs.PHI_SUB[i+1] = true
    is_free_fs.DPHI_SUB[i+1] = true
    is_continued_fs.BETA_SUB[i+1] = false
    is_continued_fs.PHI_SUB[i+1] = false
    is_continued_fs.DPHI_SUB[i+1] = false
    is_continued_fs.N_SUB[i+1] = false
  end

  is_free_fs.BETA_SUB[inN_SUB_length+1] = false
  is_free_fs.PHI_SUB[inN_SUB_length+1] = false
  is_free_fs.DPHI_SUB[inN_SUB_length+1] = true
  is_continued_fs.BETA_SUB[inN_SUB_length+1] = false
  is_continued_fs.PHI_SUB[inN_SUB_length+1] = false
  is_continued_fs.DPHI_SUB[inN_SUB_length+1] = false

  return initial_fs, is_free_fs, is_continued_fs

end

find_pert_to_double = function(inB, inL, inFROUDE, inPHI_LIM, inN, inCLUSTER_VAL)
  --[[ This starts by finding a uniform stream solution with the desired box
         width. It then continues via box height and eta0 until it hits a 'limit'
         box height which is very slightly negative, and small. Then it heads via
         fixed box height, but varying F, to a double hump solution with a
         final Froude number. ]]--

  --[[ \todo : fix up the PHI_LIM issue with F=1.10 (should be around 12-13) ]]--

  local initial_FROUDE  = 1.10
  -- local width_factor    = 0.65 / inL
  local limitA          = -0.01
  local helperETA0      = 0.05
  local display_initial = false

  io.write('\n---------------------------------------\n')
  io.write('Running script FIND_DOUBLELIM_FROM_UNIF\n')
  io.write('---------------------------------------\n')
  io.write('Begin with a near uniform-stream solution with smoothed box disturbance.\n')
  io.write('Topography parameters, B = ', inB, ', L = ', inL, ', and continue to\n')
  io.write('to a double sol. wave with A = ', limitA, '. Then find a double sol. wave solution\n')
  io.write('with F = ', inFROUDE, '.\n')
  io.flush()

  local approx_LAMBDA = math.sqrt((-5/4) + 
                                  math.sqrt((25/16) + 
                                            (15/2)*(initial_FROUDE*initial_FROUDE - 1)))

  io.write('Solved approximate equation (when F = ', initial_FROUDE, ') ; LAMBDA = ', approx_LAMBDA, ' + h.o.t. \n\n')

  local inDPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- A and D are 'magic' numbers that get us started near a uniform stream
  local initial_fs = { A = 1e-5,
                       B = inB,
                       L = inL,
                       FROUDE = initial_FROUDE,
                       LAMBDA = approx_LAMBDA,
                       PHI_C = inL/2.0,
                       U_MIN = 1.0,
                       D = 0.0005,
                       THETA_S = {},
                       BETA_SUB = { 1, inN },
                       PHI_SUB = {0, inPHI_LIM},
                       DPHI_SUB = {inDPHI_LIM/inN, inDPHI_LIM/inN},
                       N_SUB = { inN-1 },
                       CLUSTER_VAL = inCLUSTER_VAL,
                       HOMOTOPY_S = 1.0,
                       GAMMA = 25.0, -- see note??
                       ETA0 = 0
                     }

  for j = 1, inN do
    initial_fs.THETA_S[j] = 0
  end

  local _, is_free_fs, is_continued_fs = pert_to_double_tables(inN, 1)

  io.write('  Creating new surface object ... ')
  io.flush()
  local freesurface = surface.new(initial_fs,
                                  is_free_fs,
                             is_continued_fs,
           htopography_residuals.residual_hs,
     htopography_residuals.compute_output_hs)
  io.write('done.\n')

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
  io.write('         N :')
  for k,v in pairs(initial_unif.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('         A : ', initial_unif.A, ',\n')
  io.write('         B : ', initial_unif.B, ',\n')
  io.write('         D : ', initial_unif.D, ',\n')
  io.write('    LAMBDA : ', initial_unif.LAMBDA, ',\n')
  io.write('      ETA0 : ', initial_unif.ETA0, ',\n')
  io.write('     PHI_C : ', initial_unif.PHI_C, '.\n\n')

  if res ~= continuator.CONVERGED then
     print_out = true
     free_surface.residual(freesurface.unknowns, {freesurface})
    -- gtopography_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find initial guassian topography solution.')
  end

  io.write('Begin continuation with respect to A.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 100
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.6
  local started_converge = false
  local default_ds = 0.03
  local ds = default_ds
  local max_stepsize = 2.0
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
    io.write('        A : ', csurf.A, '\n')
    io.write('     ETA0 : ', csurf.ETA0, '\n')
    io.write('    PHI_C : ', csurf.PHI_C, '\n')

    if (not initial_step) and
       (csurf.A < 0) then
      ds = default_ds * 0.3 -- 0.0075
    else
      ds = default_ds
    end

    if (not initial_step) and
       math.abs(csurf.A - limitA) and
       (csurf.A < 0) and
       (csurf.ETA0 < (initial_FROUDE-1)) and
       (csurf.ETA0 > 0)
      then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
	 end
      res, msg =  freesurface:progress_approach(limitA, rough_converge_rate, freesurface)
      io.write(msg)

    else

      res, msg = freesurface:progress_interp(ds, direction, freesurface)
      io.write(msg)

    end

    if res == continuator.BIFURCATION then
      io.write('    Bifurcation detected, attempting to switch.\n')

      if DIRECTION == continuator.FORWARD then
        DIRECTION = continuator.BACKWARD
      else
        DIRECTION = continuator.FORWARD
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
 
      if csurf.A > nsurf.A then
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
    else
      freesurface, direction = htopography_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end

  local initial_help = surface.unpack_vec(freesurface.unknowns,
                                           freesurface.initial,
                                           freesurface.is_free,
                                      freesurface.is_continued,
                                         freesurface.ind_table)
 
  initial_help.A = limitA

  io.write('  Parameters for double hump solution with small A,\n')
  io.write('         N : ')
  for k,v in pairs(initial_help.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('         A : ', initial_help.A, ',\n')
  io.write('         B : ', initial_help.B, ',\n')
  io.write('         D : ', initial_help.D, ',\n')
  io.write('    LAMBDA : ', initial_help.LAMBDA, ',\n')
  io.write('      ETA0 : ', initial_help.ETA0, ',\n')
  io.write('     PHI_C : ', initial_help.PHI_C, '.\n\n')

  local is_free_help = htopography_grid.deepcopy(freesurface.is_free)
  local is_continued_help = htopography_grid.deepcopy(freesurface.is_continued)

  is_free_help.A = false
  is_free_help.FROUDE = true
  is_continued_help.A = false
  is_continued_help.FROUDE = true

  io.write('  Creating new surface object...')
  io.flush()

  freesurface = surface.new(   initial_help,
                               is_free_help,
                          is_continued_help,
          htopography_residuals.residual_hs,
    htopography_residuals.compute_output_hs)
          
  io.write('done.\n')

  io.write('  Attempting to solve with new continuation variables...')
  io.flush()
  local res = freesurface:progress(0, DIRECTION, freesurface)
  io.write('done.\n')


  if res ~= continuator.CONVERGED then
     print_out = true
     htopography_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find near unforced solitary wave solution.')
  end

  io.write('Begin continuation with respect to F.\n')

  freesurface:set_length_vars({'FROUDE', 'ETA0'})

  num_continuate = 100
  stepsize = 0.05
  direction = continuator.FORWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  default_ds = 0.01
  ds = default_ds -- 0.01 a bit too small?
  max_stepsize = 2.5
  initial_step = true

  freesurface:set_max_stepsize(max_stepsize)
  freesurface:set_stepsize(stepsize)

  res = nil
  msg = nil

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
    io.write('         A : ', csurf.A, '\n')
    io.write('    FROUDE : ', csurf.FROUDE, '\n')
    io.write('      ETA0 : ', csurf.ETA0, '\n')
    io.write('     PHI_C : ', csurf.PHI_C, '\n')

    if (not initial_step) and
       (math.abs(csurf.FROUDE - inFROUDE) < (1.2*ds)) then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
      end
      res, msg =  freesurface:progress_approach(inFROUDE, rough_converge_rate, freesurface)
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
      if inFROUDE <  csurf.FROUDE then
        lsgn = -1
      end
      if lsgn*(csurf.FROUDE) > lsgn*(nsurf.FROUDE) then
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

    else
      freesurface, direction = htopography_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end


  local initial_dbl = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)

  initial_dbl.FROUDE = inFROUDE

  local is_free_dbl = htopography_grid.deepcopy(freesurface.is_free)
  local is_continued_dbl = htopography_grid.deepcopy(freesurface.is_continued)

  is_free_dbl.A = true
  is_free_dbl.FROUDE = false
  is_continued_dbl.A = true
  is_continued_dbl.FROUDE = false

  return initial_dbl, is_free_dbl, is_continued_dbl

end
