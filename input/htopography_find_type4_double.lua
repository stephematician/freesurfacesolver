type4_double_tables = function(inN, inN_SUB_length)

  local initial_fs = { A = 0,
                       B = 0,
                       L = 0,
                       FROUDE = 0,
                       LAMBDA = 0,
                       PHI_C = 0,
                       U_MIN = 0,
                       D = 0,
                       THETA_S = {},
                       BETA_SUB = {0, 0},
                       PHI_SUB = {0, 0},
                       DPHI_SUB = {0, 0},
                       N_SUB = {0},
                       CLUSTER_VAL = 0,
                       HOMOTOPY_S = 0,
                       GAMMA = 0, -- see note??
                       ETA0 = 0
                     }

  for j = 1, inN do
    initial_fs.THETA_S[j] = 0
  end

  local is_free_fs = { A = false,
                       B = false,
                       L = false,
                       FROUDE = true,
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

  for i = 1, inN do
    is_free_fs.THETA_S[i] = true
  end
  is_free_fs.THETA_S[1] = false

  local is_continued_fs = { A = false,
                            B = false,
                            L = false,
                            FROUDE = false,
                            LAMBDA = false,
                            PHI_C = false,
                            U_MIN = true,
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

find_type4_double = function(inA, inB, inL, inU_MIN, inPHI_LIM, inN, inCLUSTER_VAL, inFROUDE)
  --[[ This starts by finding a uniform stream solution with the desired width.
       It then continues via A and ETA0 until it hits an intermediate A which
       is very slightly negative, and small. Then it heads via fixed A, but
       varying F, to a single solitary wave solution with the desired U at the
       crest then it continues with fixed U to the desired Froude number.

       Final free variables are - fixed U, varying F and A
                                - continuing F --]]

  local helperFROUDE    = 1.10
  if inFROUDE then
    helperFROUDE = inFROUDE
  end
  local helperETA0      = 1.5*(helperFROUDE - 1)
  local display_initial = false

  io.write('\n-------------------------------------\n')
  io.write(' Running script FIND_STOKES_DOUBLE \n')
  io.write('-------------------------------------\n')
  io.write('Begin with a near uniform-stream solution with smoothed box disturbance.\n')
  io.write('Topography parameters, B = ', inB, ', L = ', inL, ', and continue to\n')
  io.write('to a double wave with A = ', inA, '. Then find a  dbl. wave solution with\n')
  io.write('U_MIN = ', inU_MIN, '.\n')
  io.flush()

  local approx_LAMBDA = math.sqrt((-5/4) + 
                                  math.sqrt((25/16) + 
                                            (15/2)*(helperFROUDE*helperFROUDE - 1)))

  io.write('Solved approximate equation (when F = ', helperFROUDE, ') ; LAMBDA = ', approx_LAMBDA, ' + h.o.t. \n\n')

  local inDPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- A and D are 'magic' numbers that get us started near a uniform stream
  local initial_fs = { A = 1e-5,
                       B = inB,
                       L = inL,
                       FROUDE = helperFROUDE,
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

  local is_free_fs = { A = true,
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

  for i = 1, inN do
    is_free_fs.THETA_S[i] = true
  end
  is_free_fs.THETA_S[1] = false

  local is_continued_fs = { A = true,
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
  for i = 1, inN do
    is_continued_fs.THETA_S[i] = false
  end

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

     error('Could not find initial guassian topography solution.')
  end

  io.write('Begin continuation with respect to A.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 100
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.6
  local started_converge = false
  local default_ds = 0.05
  local ds = default_ds
  local max_stepsize = 2.5
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

--    if (not initial_step) and
--       (csurf.A < 0) then
--      ds = default_ds * 0.3 -- 0.0075
--    else
--      ds = default_ds
--    end

    if (not initial_step) and
       (math.abs(csurf.A - inA) < 1.5*ds) and
       (csurf.ETA0 > 0) and
       ((csurf.A * inA) > 0)
      then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
	 end
      res, msg =  freesurface:progress_approach(inA, rough_converge_rate, freesurface)
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
 
      if nsurf.A < csurf.A then
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
 
  initial_help.A = inA

  io.write('  Parameters for helper solitary wave solution with small A,\n')
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

  local _,is_free_help, is_continued_help = type4_double_tables(inN, #initial_help.N_SUB);

  --local is_free_help = htopography_grid.deepcopy(freesurface.is_free)
  --local is_continued_help = htopography_grid.deepcopy(freesurface.is_continued)

  --is_free_help.A = false
  --is_free_help.FROUDE = true
  --is_continued_help.A = false
  --is_continued_help.U_MIN = true

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

     error('Could not find helper solitary wave solution.')
  end

  io.write('Begin continuation with respect to U_MIN.\n')

  freesurface:set_length_vars({'U_MIN', 'ETA0'})

  num_continuate = 100
  stepsize = 0.05
  direction = continuator.BACKWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  default_ds = 0.015
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
    io.write('     U_MIN : ', csurf.U_MIN, '\n')
    io.write('    FROUDE : ', csurf.FROUDE, '\n')
    io.write('      ETA0 : ', csurf.ETA0, '\n')
    io.write('     PHI_C : ', csurf.PHI_C, '\n')

    if (not initial_step) and
       (math.abs(csurf.U_MIN - inU_MIN) < (1.2*ds)) then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
      end
      res, msg =  freesurface:progress_approach(inU_MIN, rough_converge_rate, freesurface)
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
      if (nsurf.ETA0 > csurf.ETA0) then
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
      -- freesurface, direction = htopography_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end

  freesurface, _ = htopography_grid.regrid_hs(freesurface)

  local initial_dbl = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)
 
  initial_dbl.U_MIN = inU_MIN

  local is_free_dbl = htopography_grid.deepcopy(freesurface.is_free)
  local is_continued_dbl = htopography_grid.deepcopy(freesurface.is_continued)

  return initial_dbl, is_free_dbl, is_continued_dbl

end
