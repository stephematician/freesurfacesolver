find_pert_to_double = function(inB,  inU_MIN, inPHI_LIM, inN, inCLUSTER_VAL, inA, inFROUDE)

  --[[ This starts by finding a uniform stream solution with the value of B and
       with F = 1.1. It then continues via A until it finds a double hump
       solution with A = helperA. Then it fixes A (and frees F), and traverses via U_MIN
       down to a double hump solution with the required finalU_MIN.
       After this, it traverses to a very small A = limitA by fixing U_MIN
       and free'ing A. --]]

  local initial_FROUDE  = 1.10
  local helperA         = -0.01 * inB
  if inFROUDE then
    initial_FROUDE = inFROUDE
  end
  if inA then
    helperA = inA
  end

  local display_initial = false

  io.write('\n---------------------------------------\n')
  io.write('Running script FIND_DOUBLELIM_FROM_UNIF\n')
  io.write('---------------------------------------\n')
  io.write('Begin with a near uniform-stream solution with B = ', inB, ', and F = ', initial_FROUDE, '.\n')
  io.write('Continue from this solution to a double sol. wave with A = ', helperA, '.\n')
  io.write('Then proceed from there to a double sol. wave solution with ')
  io.write('U_MIN = ', inU_MIN, '.\n\n')
  io.flush()

  local deepcopy = function(object)
    local lookup_table = {}
    local function _copy(object)
      if type(object) ~= "table" then
        return object
      elseif lookup_table[object] then
         return lookup_table[object]
      end
      local new_table = {}
      lookup_table[object] = new_table
      for index, value in pairs(object) do
        new_table[_copy(index)] = _copy(value)
      end
      return setmetatable(new_table, getmetatable(object))
    end
    return _copy(object)
  end

  local approx_LAMBDA = math.sqrt((-5/4) + 
                                  math.sqrt((25/16) + 
                                            (15/2)*(initial_FROUDE*initial_FROUDE - 1)))

  io.write('\nSolved approximate equation; LAMBDA = ', approx_LAMBDA, ' + h.o.t. \n\n')

  local inDPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- BOX_HEIGHT and Du, Dv are 'magic' numbers that get us started near a uniform stream
  local initial_surf = { A = 0.0001,
                         B = inB,
                         FROUDE = initial_FROUDE,
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
    initial_surf.V_S[j] = 0
  end

  local is_free_surf = { A = true,
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
    is_free_surf.V_S[i] = true
  end
  is_free_surf.V_S[1] = false

  --[[ The next table 'is_continued' specifies which variables are to be continued
     from. At the moment ONE and ONLY ONE value may be set to true, the rest must
     be false --]]

  local is_continued_surf = { A = true,
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
    is_continued_surf.V_S[i] = false
  end

  io.write('  Creating new surface object ... ')
  io.flush()

  local freesurface = surface.new(initial_surf,
                                  is_free_surf,
                             is_continued_surf,
                pressure_residuals.residual_hs,
          pressure_residuals.compute_output_hs)

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
  io.write('       N : ')
  for k,v in pairs(initial_unif.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       A : ', initial_unif.A, ',\n')
  io.write('       B : ', initial_unif.B, ',\n')
  io.write('       D : ', initial_unif.D, ',\n')
  io.write('       F : ', initial_unif.FROUDE, ',\n')
  io.write('  LAMBDA : ', initial_unif.LAMBDA, ',\n')
  io.write('    ETA0 : ', initial_unif.ETA0, '.\n\n')

  if res ~= continuator.CONVERGED then
     print_out_box = true
     box_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find initial near uniform solution.')
  end

  io.write('Begin continuation with respect to \'A\'.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 100
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.55
  local started_converge = false
  local ds = 0.035 -- 0.02 works ... trying 0.03
  local max_stepsize = 3.5
  local initial_step = true

  local big_ds = false

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
    io.write('          A : ', csurf.A, '\n')
    io.write('       ETA0 : ', csurf.ETA0, '\n')

    if (not initial_step) and
       (csurf.A < 0) then
      ds = 0.01 -- 0.0075
    end

    if (not initial_step) and
       (math.abs(csurf.A - helperA) < 1.5*ds) and
       (csurf.A < 0) and
       (csurf.ETA0 < (initial_FROUDE-1)) and
       (csurf.ETA0 > 0)
      then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
	 end
      res, msg =  freesurface:progress_approach(helperA, rough_converge_rate, freesurface)
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
      freesurface, direction = pressure_grid.regrid_hs(freesurface)
    end

    if inA then
      if started_converge and (i - started_converge_i) > 11 then
        io.write('  Finished.\n')
        break
      end
    else
      if started_converge and (i - started_converge_i) > 8 then
        io.write('  Finished.\n')
        break
      end
    end

  end

  local initial_dbl = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)
 
  initial_dbl.A = helperA

  io.write('  Parameters for double hump solution with A = ', helperA)
  io.write(' and F = ', initial_FROUDE, '.\n')
  io.write('       N : ')
  for k,v in pairs(initial_dbl.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       A : ', initial_dbl.A, ',\n')
  io.write('       B : ', initial_dbl.B, ',\n')
  io.write('       D : ', initial_dbl.D, ',\n')
  io.write('       F : ', initial_dbl.FROUDE, ',\n')
  io.write('  lambda : ', initial_dbl.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_dbl.ETA0, '.\n\n')

  local is_free_dbl = deepcopy(freesurface.is_free)
  local is_continued_dbl = deepcopy(freesurface.is_continued)

  is_free_dbl.A = false
  is_free_dbl.FROUDE = true
  is_continued_dbl.A = false
  is_continued_dbl.U_MIN = true

  io.write('  Creating new surface object...')
  io.flush()

  freesurface = surface.new(   initial_dbl,
                               is_free_dbl,
                          is_continued_dbl,
            pressure_residuals.residual_hs,
      pressure_residuals.compute_output_hs)
  
  io.write('done.\n')

  io.write('  Attempting to solve with new continuation variables...')
  io.flush()
  local res = freesurface:progress(0, DIRECTION, freesurface)
  io.write('done.\n')


  if res ~= continuator.CONVERGED then
     print_out = true
     pressure_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find double solitary wave solution in new continuation vars.')
  end

  io.write('Begin continuation with respect to U_MIN.\n')

  freesurface:set_length_vars({'ETA0', 'U_MIN'})

  num_continuate = 100
  stepsize = 0.05
  direction = continuator.BACKWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  ds = 0.03 -- 0.01 a bit too small?
  max_stepsize = 3.5
  initial_step = true

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
    io.write('             FROUDE : ', csurf.FROUDE, '\n')
    io.write('               ETA0 : ', csurf.ETA0, '\n')
    io.write('              U_MIN : ', csurf.U_MIN, '\n')

    if (not initial_step) and
       (math.abs(csurf.U_MIN - inU_MIN) < (2*ds)) then

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
      if (csurf.U_MIN) < (nsurf.U_MIN) then
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
      freesurface, direction = pressure_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end

  local initial_dbl_lim = surface.unpack_vec(freesurface.unknowns,
                                              freesurface.initial,
                                              freesurface.is_free,
                                         freesurface.is_continued,
                                            freesurface.ind_table)

  initial_dbl_lim.U_MIN = inU_MIN

  local is_free_dbl_lim = deepcopy(freesurface.is_free)
  local is_continued_dbl_lim = deepcopy(freesurface.is_continued)

  is_free_dbl_lim.A = true
  is_free_dbl_lim.U_MIN = false
  is_continued_dbl_lim.A = true
  is_continued_dbl_lim.U_MIN = false

  return initial_dbl_lim, is_free_dbl_lim, is_continued_dbl_lim

end
