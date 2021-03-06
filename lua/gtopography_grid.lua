require "gtopography_residuals"
require "grid_functions"

module("gtopography_grid", package.seeall)

-- Okay this has literally *one* change from gtopography_grid ... maybe I need to fix this structure.

deepcopy = function(object)
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

regrid_hs = function(insurf)

  local initial_surf = surface.unpack_vec(insurf.unknowns,
                                           insurf.initial,
                                           insurf.is_free,
                                      insurf.is_continued,
                                         insurf.ind_table)


  local is_free_surf = deepcopy(insurf.is_free)
  local is_continued_surf = deepcopy(insurf.is_continued)
  
  local oldN_SUB = deepcopy(initial_surf.N_SUB)

  local inN = #(initial_surf.THETA_S)
  local inCLUSTER_DPHI_VAL = initial_surf.CLUSTER_VAL
  local inPHI_LIM = initial_surf.PHI_SUB[#initial_surf.PHI_SUB]

  local ibeta_s, itheta_s_mid, itau_s_mid,
                     ix_s_mid,     ix_lhs,
                   ieta_s_mid,   ieta_lhs = 
                            gtopography_residuals.compute_bttxe_hs(initial_surf)

  local iphi_s  = integrals_sbox.phi(ibeta_s,
                       initial_surf.BETA_SUB,
                        initial_surf.PHI_SUB,
                       initial_surf.DPHI_SUB,
                     initial_surf.HOMOTOPY_S)

   initial_surf.THETA_S,
  initial_surf.BETA_SUB,
     initial_surf.N_SUB,
   initial_surf.PHI_SUB,
  initial_surf.DPHI_SUB = grid_functions.convert_grid_hs(iphi_s,
                                           initial_surf.THETA_S,
                                                     itau_s_mid,
                                                            inN,
                                             inCLUSTER_DPHI_VAL,
                                        initial_surf.HOMOTOPY_S,
                                             integrals_sbox.phi)

  local newN_SUB = initial_surf.N_SUB

  local need_regrid = false
  if #newN_SUB == #oldN_SUB then
    for k = 1, #newN_SUB do
      if not (newN_SUB[k] == oldN_SUB[k]) then
         need_regrid = true
      end
    end
  else
    need_regrid = true
  end

  if need_regrid then
    io.write('  Performing regrid ... ')
    io.flush()

    local numN_SUB = #(initial_surf.N_SUB)

    is_free_surf.PHI_SUB = { false }
    is_free_surf.BETA_SUB = { false }
    is_free_surf.DPHI_SUB = { true }
    is_free_surf.N_SUB = { false }
    is_continued_surf.PHI_SUB = { false }
    is_continued_surf.BETA_SUB = { false }
    is_continued_surf.DPHI_SUB = { false }
    is_continued_surf.N_SUB = { false }

    for i = 1, #newN_SUB-1 do
      is_free_surf.BETA_SUB[i+1] = false
      is_free_surf.N_SUB[i+1] = false
      is_free_surf.PHI_SUB[i+1] = true
      is_free_surf.DPHI_SUB[i+1] = true
      is_continued_surf.BETA_SUB[i+1] = false
      is_continued_surf.PHI_SUB[i+1] = false
      is_continued_surf.DPHI_SUB[i+1] = false
      is_continued_surf.N_SUB[i+1] = false
    end

    is_free_surf.BETA_SUB[#newN_SUB+1] = false
    is_free_surf.PHI_SUB[#newN_SUB+1] = false
    is_free_surf.DPHI_SUB[#newN_SUB+1] = true
    is_continued_surf.BETA_SUB[#newN_SUB+1] = false
    is_continued_surf.PHI_SUB[#newN_SUB+1] = false
    is_continued_surf.DPHI_SUB[#newN_SUB+1] = false

    local newsurface = surface.new(initial_surf,
                                   is_free_surf,
                              is_continued_surf,
                                insurf.residual,
                          insurf.compute_output)

    local c_direction = insurf:get_direction()

    local res = newsurface:progress(0, c_direction, newsurface)

    -- This doesn't account for bifurcation :(
    assert(res == continuator.CONVERGED, '... Could not regrid, failed to' ..
           ' converge with new grid.')

    local csurf = surface.unpack_vec(newsurface.unknowns,
                                      newsurface.initial,
                                      newsurface.is_free,
                                  newsurface.is_continued,
                                     newsurface.ind_table)

    newsurface:set_max_stepsize(insurf:get_max_stepsize())
    newsurface:set_prev_stepsize(insurf:get_prev_stepsize())
    newsurface:set_prev_length_vars(insurf:get_prev_length_vars())
    newsurface:set_length_vars(insurf:get_length_vars())
    newsurface:set_prev_finalval(insurf:get_prev_finalval())

    res = newsurface:progress(insurf:get_stepsize() * 0.01,
                                               c_direction,
                                                newsurface)
    -- This doesn't account for bifurcation :(
    assert(res == continuator.CONVERGED, '... Could not regrid, failed to' ..
           ' converge with new grid.')
    io.write('done.\n')
    io.flush()

    local nsurf = surface.unpack_vec(newsurface.unknowns,
                                      newsurface.initial,
                                      newsurface.is_free,
                                  newsurface.is_continued,
                                     newsurface.ind_table)

    local switched = newsurface:test_switch(nsurf, csurf, insurf:get_prev_length_vars())

    if switched then
      io.write('  Regrid caused direction change.\n')
      io.flush()
      if c_direction == continuator.FORWARD then
        c_direction = continuator.BACKWARD
      else
        c_direction = continuator.FORWARD
      end

      res = newsurface:progress(insurf:get_stepsize() * 0.01,
                                      c_direction,
                                      newsurface)
      -- This doesn't account for bifurcation :(
      assert(res == continuator.CONVERGED, '... Could not regrid, failed to' ..
           ' converge with new grid.')
    end

    return newsurface, c_direction -- deepcopy(insurf.prev_length_vars)
  else
    return insurf, insurf:get_direction()
  end

end
