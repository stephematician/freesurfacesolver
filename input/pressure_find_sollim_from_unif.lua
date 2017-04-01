find_pert_to_sol = function(inB, finalFROUDE, inU_MIN, inPHI_LIM, inN, inCLUSTER_VAL)
  --[[ This starts by finding a uniform stream solution with the desired box
       width. It then continues via box height and eta0 until it hits a 'limit'
       box height which is very slightly negative, and small. Then it heads via
       fixed box height, but varying F, to a double hump solution with a
       final Froude number. --]]

   --[[ \todo : fix up the PHI_LIM issue with F=1.10 (should be around 12-13) --]]

  local initial_FROUDE  = 1.10
  local display_initial = false

  io.write('\n---------------------------------------\n')
  io.write('Running script FIND_SOLLIM_FROM_UNIF\n')
  io.write('---------------------------------------\n')
  io.write('Begin with a near uniform-stream solution with B = ', inB, ', and \n')
  io.write('F = ', initial_FROUDE, ', then continue solutions to a solitary wave with\n')
  io.write('U_MIN = ', inU_MIN, '. Then find a solitary wave solution with ')
  io.write('F = ', finalFROUDE, '.\n')
  io.flush()

  helperA = -0.0001 * inB

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

  io.write('Solved approximate equation for lambda : ', approx_LAMBDA, '\n\n')

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
     be false ]]--

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
  io.write('             n : ')
  for k,v in pairs(initial_unif.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       A : ', initial_unif.A, ',\n')
  io.write('       B : ', initial_unif.B, ',\n')
  io.write('       D : ', initial_unif.D, ',\n')
  io.write('  lambda : ', initial_unif.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_unif.ETA0, '.\n\n')

  if res ~= continuator.CONVERGED then
     print_out_box = true
     box_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find initial near uniform solution.')
  end

  io.write('Begin continuation with respect to \'a\'.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 100
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.55
  local started_converge = false
  local ds = 0.035 -- 0.02 works ... trying 0.03
  local max_stepsize = 2.5
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
    io.write('          a : ', csurf.A, '\n')
    io.write('       eta0 : ', csurf.ETA0, '\n')

    --[[if (not initial_step) and
       (csurf.A < 0) then
      ds = 0.0075
    end --]]

    if (not initial_step) and
       math.abs(csurf.A - helperA) < ds and
       (csurf.ETA0 > (initial_FROUDE-1))
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

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end
  local initial_unf = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)
 
  initial_unf.A = helperA

  io.write('  Parameters for near-unforced single solitary wave solution,\n')
  io.write('             n : ')
  for k,v in pairs(initial_unf.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       a : ', initial_unf.A, ',\n')
  io.write('       b : ', initial_unf.B, ',\n')
  io.write('       D : ', initial_unf.D, ',\n')
  io.write('  lambda : ', initial_unf.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_unf.ETA0, '.\n\n')

  local is_free_unf = deepcopy(freesurface.is_free)
  local is_continued_unf = deepcopy(freesurface.is_continued)

  is_free_unf.A = false
  is_free_unf.FROUDE = true
  is_continued_unf.A = false
  is_continued_unf.U_MIN = true

  io.write('  Creating new surface object...')
  io.flush()

  freesurface = surface.new(   initial_unf,
                               is_free_unf,
                          is_continued_unf,
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

     error('Could not find near unforced solitary wave solution.')
  end

  io.write('Begin continuation with respect to F.\n')

  freesurface:set_length_vars({'FROUDE', 'U_MIN'})

  num_continuate = 100
  stepsize = 0.05
  direction = continuator.FORWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  ds = 0.04
  max_stepsize = 2.5
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
    io.write('             froude : ', csurf.FROUDE, '\n')
    io.write('               umin : ', csurf.U_MIN, '\n')

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
      if csurf.U_MIN < nsurf.U_MIN then
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

  local initial_sollim = surface.unpack_vec(freesurface.unknowns,
                                             freesurface.initial,
                                             freesurface.is_free,
                                        freesurface.is_continued,
                                           freesurface.ind_table)

  initial_sollim.U_MIN = inU_MIN

  io.write('  Parameters for near-unforced single solitary wave solution,\n')
  io.write('             n : ')
  for k,v in pairs(initial_sollim.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       a : ', initial_sollim.A, ',\n')
  io.write('       b : ', initial_sollim.B, ',\n')
  io.write('       D : ', initial_sollim.D, ',\n')
  io.write('  lambda : ', initial_sollim.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_sollim.ETA0, '.\n\n')

  local is_free_sollim = deepcopy(freesurface.is_free)
  local is_continued_sollim = deepcopy(freesurface.is_continued)

  is_free_sollim.A = true
  is_free_sollim.U_MIN = false
  is_continued_sollim.FROUDE = true
  is_continued_sollim.U_MIN = false

  io.write('  Creating new surface object...')
  io.flush()

  freesurface = surface.new(initial_sollim,
                            is_free_sollim,
                       is_continued_sollim,
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

     error('Could not find near unforced solitary wave solution.')
  end

  io.write('Begin continuation with respect to F.\n')

  freesurface:set_length_vars({'FROUDE', 'ETA0'})

  num_continuate = 100
  stepsize = 0.01
  direction = continuator.FORWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  ds = 0.01
  max_stepsize = 2.5
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
    io.write('             froude : ', csurf.FROUDE, '\n')
    io.write('               eta0 : ', csurf.ETA0, '\n')

    if (not initial_step) and
       (math.abs(csurf.FROUDE - finalFROUDE) < (1.5*ds)) then

      if not started_converge then
        started_converge = true
        started_converge_i = i
        io.write('  Switching step size selection routine.\n')
      end
      res, msg =  freesurface:progress_approach(finalFROUDE, rough_converge_rate, freesurface)
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
      if finalFROUDE <  csurf.FROUDE then
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
      freesurface, direction = pressure_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 10 then
      io.write('  Finished.\n')
      break
    end

  end

  local initial_solpert = surface.unpack_vec(freesurface.unknowns,
                                              freesurface.initial,
                                              freesurface.is_free,
                                         freesurface.is_continued,
                                            freesurface.ind_table)

  initial_solpert.U_MIN = inU_MIN

  io.write('  Parameters for near-unforced single solitary wave solution,\n')
  io.write('             n : ')
  for k,v in pairs(initial_sollim.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       A : ', initial_solpert.A, ',\n')
  io.write('       D : ', initial_solpert.D, ',\n')
  io.write('    eta0 : ', initial_solpert.ETA0, '.\n')
  io.write('       F : ', initial_solpert.FROUDE, '.\n')
  io.write('    u(0) : ', initial_solpert.U_MIN, '.\n\n')

  local is_free_solpert = deepcopy(freesurface.is_free)
  local is_continued_solpert = deepcopy(freesurface.is_continued)

  is_free_solpert.FROUDE = false
  is_free_solpert.U_MIN = true
  is_continued_solpert.FROUDE = false
  is_continued_solpert.A = true

  return initial_solpert, is_free_solpert, is_continued_solpert

end
