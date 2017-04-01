-- require "surface"
require "integrals_sbox_homotopy"
require "grid_functions"

module("htopography_load", package.seeall)

load_last = function(file_name, init_tables_fn, hi_res_n)

  local surf_file = assert(io.open(file_name, "r"))

  -- First step, find % 00000xx etc
  surf_file:seek("set", string.len("%% "))
  local n = surf_file:read("*n")
  
  io.write("n = " .. n .. "\n")

  -- Second step, parse every line that is of the form
  -- o.*{xx} = [ .... ];
  -- o.*{xx} = yyyy;

  local input_fs = {}

  local line = surf_file:read("*line")
  while (not (line == nil)) do
    local a, b, var_str = string.find(line, 'o%.(.+){' .. n .. '}')
    if a then
      var_str = string.upper(var_str)
      local d = b
      local count = 1
      local val = nil
      local _, c = string.find(line, '%[', d+1)
      if not (c == nil) then
        -- we have a table
        input_fs[var_str] = {}
        _, d, val = string.find(line, '(%-?%d+%.?%d*e?-?%d*)', c+1)
        while (not (d == nil)) do
          (input_fs[var_str])[count] = tonumber(val)
          _, d, val = string.find(line, '(%-?%d+%.?%d*e?-?%d*)', d+1)
          count = count + 1
        end
      else
        -- should be a single number
        _, d, val = string.find(line, '(%-?%d+%.?%d*e?-?%d*)', d+1)
        input_fs[var_str] = tonumber(val)
      end
    end
    io.flush()
    line = surf_file:read("*line")
  end

  surf_file:close()

  -- Third step, fill out a real free surface table with the values from
  -- temp_fs
  -- this should work for everything but THETA_S, which needs to
  -- be interpolated onto the correct grid and then solved for.

  -- Should really test a whole bunch of assertions here. These will do for now
  assert(not (input_fs.THETA_S_MID == nil), 'Could not load previous theta_s_mid data')
  assert(not (input_fs.PHI_S_MID == nil), 'Could not load previous theta_s_mid data')
  assert(not (input_fs.N_SUB == nil), 'Could not load previous n_sub data')
  assert(not (input_fs.BETA_SUB == nil), 'Could not load previous beta_sub data')
  assert(not (input_fs.PHI_SUB == nil), 'Could not load previous phi_sub data')
  assert(not (input_fs.DPHI_SUB == nil), 'Could not load previous dphi_sub data')
  assert(not (input_fs.HOMOTOPY_S == nil), 'Could not load previous homotopy_s data')

  local inN = #(input_fs.THETA_S_MID) + 1

  local initial_fs, _, _ = init_tables_fn(inN, #(input_fs.N_SUB))

  local ibeta_s = grid_functions.calculate_beta(input_fs.N_SUB)

  local iphi_s  = integrals_sbox.phi(ibeta_s,
                           input_fs.BETA_SUB,
                            input_fs.PHI_SUB,
                           input_fs.DPHI_SUB,
                         input_fs.HOMOTOPY_S)


  for k, v in pairs(initial_fs) do
    if type(v) == "number" then
      assert(not (input_fs[k] == nil), "Could not find corresponding old data for " .. k)
      initial_fs[k] = input_fs[k]
    else
      assert(type(v) == "table", "output data not in correct format")
      if not (k == "THETA_S") then
        assert(not (input_fs[k] == nil), "Could not find corresponding old data for " .. k)
        initial_fs[k] = htopography_grid.deepcopy(input_fs[k])
      end
    end
  end

  initial_fs.THETA_S = grid_functions.interp_cubic(input_fs.PHI_S_MID, input_fs.THETA_S_MID, iphi_s)

  initial_fs.THETA_S[1] = 0 -- given

   initial_fs.THETA_S,
  initial_fs.BETA_SUB,
     initial_fs.N_SUB,
   initial_fs.PHI_SUB,
  initial_fs.DPHI_SUB = grid_functions.convert_grid_hs(iphi_s,
                                           initial_fs.THETA_S,
                                           input_fs.TAU_S_MID,
                                                     hi_res_n,
                                       initial_fs.CLUSTER_VAL,
                                        initial_fs.HOMOTOPY_S,
                                           integrals_sbox.phi)

  local _, is_free_fs, is_continued_fs = init_tables_fn(hi_res_n, #(initial_fs.N_SUB))

  io.write('  Parameters for solution from old data,\n')
  io.write('       N :')
  for k,v in pairs(initial_fs.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('       A : ', initial_fs.A, ',\n')
  io.write('       B : ', initial_fs.B, ',\n')
  io.write('       D : ', initial_fs.D, ',\n')
  io.write('  LAMBDA : ', initial_fs.LAMBDA, ',\n')
  io.write('    ETA0 : ', initial_fs.ETA0, ',\n')
  io.write('   U_MIN : ', initial_fs.U_MIN, ',\n')
  io.write('   PHI_C : ', initial_fs.PHI_C, '.\n\n')

  io.write('  Creating new surface object from old data ... ')
  io.flush()
  local freesurface = surface.new(initial_fs,
                                  is_free_fs,
                             is_continued_fs,
           htopography_residuals.residual_hs,
     htopography_residuals.compute_output_hs)
  io.write('done.\n')

  local DIRECTION = continuator.BACKWARD

  io.write('  Attempting to sort out old grid ... ')
  io.flush()
  local res = freesurface:progress(0, DIRECTION, freesurface)
  io.write('done.\n')
  io.flush()

  loaded_fs = surface.unpack_vec(freesurface.unknowns,
                                  freesurface.initial,
                                  freesurface.is_free,
                             freesurface.is_continued,
                                freesurface.ind_table)

  io.write('  Parameters for nearest solution to old data,\n')
  io.write('       N :')
  for k,v in pairs(loaded_fs.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('       A : ', loaded_fs.A, ',\n')
  io.write('       B : ', loaded_fs.B, ',\n')
  io.write('       D : ', loaded_fs.D, ',\n')
  io.write('  LAMBDA : ', loaded_fs.LAMBDA, ',\n')
  io.write('    ETA0 : ', loaded_fs.ETA0, ',\n')
  io.write('   U_MIN : ', loaded_fs.U_MIN, ',\n')
  io.write('   PHI_C : ', loaded_fs.PHI_C, '.\n\n')

  return loaded_fs, is_free_fs, is_continued_fs

end
