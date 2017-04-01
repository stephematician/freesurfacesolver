require "gtopography_residuals"
require "sbox_residuals"

module("sbox_grid_functions", package.seeall)

--[[ Interpolate theta using cubic spline interpolation.

    arguments:
    phi   - old grid coordinates.
    theta - value of theta on old grid.
    nphi  - new grid coordinates.

    returns:
    ntheta - theta interpolated onto new grid coordinates.

    requires phi and nphi are monotonically increasing, and
    that #theta == #phi.

--]]
interp_cubic = function(phi, theta, nphi)
  local b = 0
  local ntheta = {}

  assert(#phi == #theta, "Function interp_cubic\n" ..
                         "Failed pre-condition #phi == #theta\n" ..
                         "Exit\n\n")
  for k = 2, #phi do
     assert(phi[k] > phi[k-1], "Function interp_cubic\n" ..
                        "Failed pre-condition phi is monotonic\n" ..
                        "Exit\n\n")
  end

  for k = 2, #nphi do
     assert(nphi[k] > nphi[k-1], "Function interp_cubic\n" ..
                        "Failed pre-condition nphi is monotonic\n" ..
                        "Exit\n\n")
  end

  for a = 1, #nphi do
    -- Spline interpolation code
    if not (b == #theta) then
      local cdist = phi[b+1] - nphi[a]
      while cdist < 0 do
        b = b + 1
        if (b == #theta) then break end
        cdist = phi[b+1] - nphi[a]
      end
    end

    local x, x0, x1, x2, x3
    local f0, f1, f2, f3
    if (b > 1) and (b < #theta-1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b-1], phi[b], phi[b+1], phi[b+2]
      f0, f1, f2, f3 = theta[b-1], theta[b], theta[b+1], theta[b+2]
    elseif (b == 1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b], phi[b+1], phi[b+2], phi[b+3]
      f0, f1, f2, f3 = theta[b], theta[b+1], theta[b+2], theta[b+3]
    elseif (b == #theta-1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b-2], phi[b-1], phi[b], phi[b+1]
      f0, f1, f2, f3 = theta[b-2], theta[b-1], theta[b], theta[b+1]
    end

    if not ((b == #theta) or (b == 0)) then
      ntheta[a] = f0*((x - x1)*(x - x2)*(x - x3))/((x0 - x1)*(x0 - x2)*(x0 - x3)) +
                  f1*((x - x0)*(x - x2)*(x - x3))/((x1 - x0)*(x1 - x2)*(x1 - x3)) +
                  f2*((x - x0)*(x - x1)*(x - x3))/((x2 - x0)*(x2 - x1)*(x2 - x3)) +
                  f3*((x - x0)*(x - x1)*(x - x2))/((x3 - x0)*(x3 - x1)*(x3 - x2))
    elseif b == 0 then
      ntheta[a] = theta[1]
    else
      ntheta[a] = theta[#theta]
    end

  end

  assert(#ntheta == #nphi, "Function interp_cubic\n" ..
                           "Failed post-condition #nphi == #ntheta\n" ..
                           "Exit\n\n") 

  return ntheta
end

--[[ Find peaks in free surface via linear interpolation

     finds zeros of theta by linear interpolation, and records location
     in terms of 'index' i into phi, and the value of phi(i)

     arguments:
     phi   - grid coordinates.
     theta - value of theta on old grid.

     returns:
     phi_peaks - value of phi at peaks
     i_peaks   - 'index' into grid for peak (non-integer)

     ensures that phi_peaks(i) lies in the correct interval of phi given i_peaks(i)
]]--
find_peaks = function(phi, theta)
  local phi_peaks = {}
  local i_peaks = {}

  assert(#phi == #theta, "Function find_peaks\n" ..
                         "Failed pre-condition #phi == #theta\n" ..
                         "Exit\n\n")
  --[[
  for kk, vv in pairs(phi) do
    print(kk, vv)
  end
  for kk, vv in pairs(theta) do
    print(kk, vv)
  end
  ]]

  -- Experimental stuff - need to know how 'peaky' a peak is.
 --[[ max_dtheta = 0
  for i = 2, (#theta) do
    if math.abs(theta[i] - theta[i-1]) > max_dtheta then
      max_dtheta = math.abs(theta[i] - theta[i-1])
    end
  end
  peak_cutoff = max_dtheta * 0.1]]
  -- /Experimental stuff

  for i = 2, (#theta) do
     if theta[i] * theta[i-1] < 0 and theta[i-1] > 0 then
      -- We have found a 'peak' between i-1 and i, now we
      -- determine a more 'precise' location of the peak via linear interpolation.
      local t = theta[i-1] / (theta[i-1] - theta[i])

      -- Experimental stuff
     -- if math.abs(theta[i-1] - theta[i]) < peak_cutoff then
      -- /Experimental stuff

      phi_peaks[#phi_peaks+1] = (1-t)*phi[i-1] + t*phi[i]
      i_peaks[#i_peaks+1]     = (i-1)+t

     -- Experimental stuff
    --  end
      -- /Experimental stuff

    end
  end
  
  print(' --------- peaks summary ------------')
  print('phi_peaks')
  for k,v in pairs(phi_peaks) do
    print(k,v)
  end
  print('i_peaks')
  for k,v in pairs(i_peaks) do
    print(k,v)
  end
  print(' --------- end print summary --------')
  
  assert(#phi_peaks == #i_peaks, "Function find_peaks\n" ..
                                 "Failed post-condition #phi_peaks == #i_peaks\n" ..
                                 "Exit\n\n") 
  for k = 1, #phi_peaks do
     assert(phi_peaks[k] < phi[math.ceil(i_peaks[k])], "Function find_peaks\n" ..
            "Failed post-condition phi_peaks(i) < phi(ceil(i_peaks(i)))\n" ..
            "for k = " .. tonumber(k) .. "." ..
            "Exit\n\n") 
     assert(phi_peaks[k] > phi[math.floor(i_peaks[k])], "Function find_peaks\n" ..
            "Failed post-condition phi_peaks(i) > phi(floor(i_peaks(i)))\n" ..
            "for k = " .. tonumber(k) .. "." ..
            "Exit\n\n") 
  end

  return phi_peaks, i_peaks
end


--[[ Converts grid from having no peak-data to having peak data. ]]
convert_grid_hs = function(phi, theta, tau_mid, n, cluster_coeff, homotopy_s)

  local compute_bsub = function(phi_crit_a, phi_crit_b, u_a, u_b)
     return (phi_crit_b - phi_crit_a) / math.pow(u_a + u_b,1) 
  end

  local calculate_u = function(theta, tau_mid, k)
    local tau_k = 0
    if k == 1 then
      tau_k = (3*tau_mid[k] - tau_mid[k+1]) / 2
    elseif k == #theta then
      tau_k = (3*tau_mid[k-1] - tau_mid[k-2]) / 2
    else
      tau_k = (tau_mid[k-1] + tau_mid[k]) / 2
    end
    
    return (math.cos(theta[k]) * math.exp(tau_k))
  end

  local phi_peaks, i_peaks = find_peaks(phi, theta)

  local phi_crit = { 0 }
  for k, phi_peaks_k in pairs(phi_peaks) do
    phi_crit[k+1] = phi_peaks_k
  end
  phi_crit[#phi_crit+1] = phi[#phi]
  
  local i_crit = { 1 }
  for k, i_peaks_k in pairs(i_peaks) do
    i_crit[k+1] = i_peaks_k
  end
  i_crit[#i_crit+1] = #theta

  --local phi_lim = phi[#phi]

  local weighted_bsub_total = 0
  for k = 1, #phi_crit-1 do
    local uka = calculate_u(theta, tau_mid, math.floor(i_crit[k]))
    local ukb = calculate_u(theta, tau_mid, math.floor(i_crit[k+1]))

    local weighted_bsub_k = compute_bsub(phi_crit[k], phi_crit[k+1], uka, ukb)
    --print('weighted_k', weighted_bsub_k)
    weighted_bsub_total = weighted_bsub_total + weighted_bsub_k
  end
  --print('weighted_total', weighted_bsub_total)

  --[[
  Initial entries are related to the centre of the wave phi = 0
  ]]

  local u0a = calculate_u(theta, tau_mid, i_crit[1])
  u0a = math.max(math.min(u0a,1.0),0.0)

  local nn_sub = { }
  local nbeta_sub = { 1 }
  local nphi_sub = { phi_crit[1] }
  local ndphi_sub = { (cluster_coeff * math.pow(1-u0a,2) + (phi_crit[#phi_crit] * math.pow(u0a,2)))/n }

  --[[
    Next entries are determined by values at peaks
  ]]
  local cu_bsub_k = 0
  local cu_nsub_k = 0
  for k = 2, #phi_crit-1 do
    nphi_sub[k] = phi_crit[k]

    local uja = calculate_u(theta, tau_mid, math.floor(i_crit[k-1]))
    local ujb = calculate_u(theta, tau_mid, math.floor(i_crit[k]))
    uja = math.max(math.min(uja,1.0),0.0)
    ujb = math.max(math.min(ujb,1.0),0.0)
    ndphi_sub[k] = (cluster_coeff * math.pow(1-ujb,2) + (phi_crit[#phi_crit] * math.pow(ujb,2))) / n

    local del_bsub_k = compute_bsub(phi_crit[k-1], phi_crit[k], uja, ujb)
    
    cu_bsub_k = cu_bsub_k + del_bsub_k
    -- math.round == math.floor(0.5 + )
    nn_sub[k-1] = math.floor(0.5 + (cu_bsub_k * (n-1) / weighted_bsub_total)) - cu_nsub_k
    cu_nsub_k = cu_nsub_k + nn_sub[k-1]
    nbeta_sub[k] = nbeta_sub[k-1] + nn_sub[k-1]
  end

  nphi_sub[#nphi_sub+1]   = phi_crit[#phi_crit]

  local uNa = calculate_u(theta, tau_mid, math.floor(i_crit[#i_crit-1]))
  local uNb = calculate_u(theta, tau_mid, math.floor(i_crit[#i_crit]))
  uNa = math.max(math.min(uNa,1.0),0.0)
  uNb = math.max(math.min(uNb,1.0),0.0)

  ndphi_sub[#ndphi_sub+1] = (cluster_coeff * math.pow(1-uNb,2) +  (phi_crit[#phi_crit] * math.pow(uNb,2))) / n
  local del_bsub_k = compute_bsub(phi_crit[#phi_crit-1], phi_crit[#phi_crit], uNa, uNb)

  cu_bsub_k = cu_bsub_k + del_bsub_k
  -- math.round == math.floor(0.5 + )
  if #phi_crit > 2 then
    nn_sub[#nn_sub+1]     = math.floor(0.5 + (cu_bsub_k * (n-1) / weighted_bsub_total)) - cu_nsub_k
  else
    nn_sub = {n-1}
  end
  cu_nsub_k = cu_nsub_k + nn_sub[#phi_crit-1]

  nbeta_sub[#nbeta_sub+1] = nbeta_sub[#nbeta_sub] + nn_sub[#nn_sub]

  local sum_nn_sub = 0
  for i = 1, #nn_sub do
    sum_nn_sub = sum_nn_sub + nn_sub[i]
  end

  assert(sum_nn_sub == (n-1), 'Failed assertion that sum(nn_sub) = n-1')
  assert(nbeta_sub[#nbeta_sub] == n, 'Failed assertion that nbeta_sub[end] = n')

  local beta   = gtopography_residuals.calculate_beta(nn_sub)
  -- yet to test this :
  local nphi   = integrals_sbox.phi(beta, nbeta_sub, nphi_sub, ndphi_sub, homotopy_s)

  local ntheta = interp_cubic(phi, theta, nphi)

  return ntheta, nbeta_sub, nn_sub, nphi_sub, ndphi_sub

end

--[[ \todo need to update this to now be 'callable' at any time in computation, and only
   regrids if necessary ]]--

regrid_hs = function(insurf)

  function deepcopy(object)
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
  
  --[[local beta   = gtopography_residuals.calculate_beta(nn_sub)]]--
  
  local ibeta_s, itheta_s_mid, itau_s_mid,
                     ix_s_mid,     ix_lhs,
                   ieta_s_mid,   ieta_lhs = box_residuals.compute_bttxe_hs(initial_surf)

  local iphi_s  = integrals_homotopy.phi(ibeta_s,
                           initial_surf.BETA_SUB,
                            initial_surf.PHI_SUB,
                           initial_surf.DPHI_SUB,
                         initial_surf.HOMOTOPY_S)

   initial_surf.THETA_S,
  initial_surf.BETA_SUB,
     initial_surf.N_SUB,
   initial_surf.PHI_SUB,
  initial_surf.DPHI_SUB = convert_grid_hs(iphi_s,
                            initial_surf.THETA_S,
                                      itau_s_mid,
                                             inN,
                              inCLUSTER_DPHI_VAL,
                         initial_surf.HOMOTOPY_S)

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

    local res = newsurface:progress(0, continuator.FORWARD, newsurface)

    newsurface:set_stepsize(insurf:get_stepsize() * 0.01)
    newsurface:set_max_stepsize(insurf:get_max_stepsize())

    newsurface:set_length_vars(insurf:get_length_vars())

    if res ~= continuator.CONVERGED then
      error('... Could not regrid')
    end

    return true, newsurface, deepcopy(insurf.prev_length_vars)
  else
    return false, insurf, nil
  end

end
