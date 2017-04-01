find_unforced_sol = function(inPHI_LIM, inB, inN, inCLUSTER_VAL, finalFROUDE)

  --[[ This starts by finding a uniform stream solution with the desired box
       width. It then continues via box height and eta0 until it hits an
       unforced solitary wave with F=1.1. Then it heads via
       fixed A, but varying F, to an unforced solitary wave with
       final Froude number. ]]--

  --[[ \todo : fix up the PHI_LIM issue with F=1.10 (should be around 12-13) ]]--

  local initial_FROUDE  = 1.10
  local width_factor    = 0.65 / inB
  local limitA          = 1e-9
  local display_initial = false

  io.write('\n-----------------------------------------\n')
  io.write('Running script FIND_UNFORCED_SOL\n')
  io.write('-----------------------------------------\n')
  io.write('Begin with a near uniform-stream solution with b = ', inB, ', and \n')
  io.write('Froude number = ', initial_FROUDE, ', the continue solutions to an unforced sol. wave with\n')
  io.write('a = ', limitA, '. Then find a unforced sol. wave solution with\n')
  io.write('Froude number = ', finalFROUDE, '.\n')
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

  io.write('Solved approximate equation for lambda : ', approx_LAMBDA, '\n\n')

  local inDPHI_LIM = 3*math.pow(PHI_LIM, 2/3)

  -- BOX_HEIGHT and Du, Dv are 'magic' numbers that get us started near a uniform stream
  local initial_surf = { A = 0.0001,
                         B = inB,
                         FROUDE = initial_FROUDE,
                         LAMBDA = approx_LAMBDA,
                         U_MIN = 1.0,
                         Du = 0.0005,
                         Dv = 0.0005,
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
                         Du = true,
                         Dv = true,
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
                              Du = false,
                              Dv = false,
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
  io.write('       a : ', initial_unif.A, ',\n')
  io.write('       b : ', initial_unif.B, ',\n')
  io.write('      Du : ', initial_unif.Du, ',\n')
  io.write('      Dv : ', initial_unif.Dv, ',\n')
  io.write('  lambda : ', initial_unif.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_unif.ETA0, '.\n\n')

  if res ~= continuator.CONVERGED then
     print_out_box = true
     box_residuals.residual_hs(freesurface.unknowns, {freesurface})

     error('Could not find initial near uniform solution.')
  end

  io.write('Begin continuation with respect to \'a\'.\n')

  freesurface:set_length_vars({'A', 'ETA0'})

  local num_continuate = 150
  local stepsize = 0.05
  local direction = continuator.FORWARD

  local rough_converge_rate = 0.55
  local started_converge = false
  local ds = 0.035 -- 0.02 works ... trying 0.03
  local max_stepsize = 10
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

    if (not initial_step) and
       (csurf.A < 0) then
      ds = 0.0075
    end

    if (not initial_step) and
       (math.abs(csurf.A - limitA) < 1.5*ds) and
       (csurf.ETA0 > (initial_FROUDE-1)) then

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
      freesurface, direction = pressure_grid.regrid_hs(freesurface)
    end

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end

  local initial_sol = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)
 
  initial_sol.A = limitA

  io.write('  Parameters for unforced solution,\n')
  io.write('             n : ')
  for k,v in pairs(initial_sol.N_SUB) do
     io.write(tostring(v) .. ' ')
  end
  io.write(',\n')
  io.write('       a : ', initial_sol.A, ',\n')
  io.write('       b : ', initial_sol.B, ',\n')
  io.write('      Du : ', initial_sol.Du, ',\n')
  io.write('      Dv : ', initial_sol.Dv, ',\n')
  io.write('  lambda : ', initial_sol.LAMBDA, ',\n')
  io.write('    eta0 : ', initial_sol.ETA0, '.\n\n')

  local is_free_sol = deepcopy(freesurface.is_free)
  local is_continued_sol = deepcopy(freesurface.is_continued)

  is_free_sol.A = false
  is_free_sol.FROUDE = true
  is_continued_sol.A = false
  is_continued_sol.FROUDE = true

  io.write('  Creating new surface object...')
  io.flush()

  freesurface = surface.new(initial_sol,
                            is_free_sol,
                       is_continued_sol,
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
  stepsize = 0.05
  direction = continuator.FORWARD

  rough_converge_rate = 0.55
  started_converge = false
  started_converge_i = 0
  ds = 0.01 -- 0.01 a bit too small?
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

    if csurf.U_MIN < 0.15 then
      ds = 0.0025
    elseif csurf.U_MIN < 0.40 then
      ds = 0.0050
    end

    if (not initial_step) and
       (math.abs(csurf.FROUDE - finalFROUDE) < (1.2*ds)) then

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

    if started_converge and (i - started_converge_i) > 9 then
      io.write('  Finished.\n')
      break
    end

  end


  local initial_unfsol = surface.unpack_vec(freesurface.unknowns,
                                             freesurface.initial,
                                             freesurface.is_free,
                                        freesurface.is_continued,
                                           freesurface.ind_table)

  initial_unfsol.FROUDE = finalFROUDE

  local is_free_unfsol = deepcopy(freesurface.is_free)
  local is_continued_unfsol = deepcopy(freesurface.is_continued)

  is_free_unfsol.A = false
  is_free_unfsol.FROUDE = true
  is_continued_unfsol.A = false
  is_continued_unfsol.FROUDE = false
  is_continued_unfsol.U_MIN = true

  return initial_unfsol, is_free_unfsol, is_continued_unfsol

end
