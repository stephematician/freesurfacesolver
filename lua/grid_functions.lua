local grid_functions = {}

--[[ Calculate value of beta for grid

     beta has been reduced to a cumulative count of how many grid points
     there are. this seems redundant as it must act like an index, but
     i'd have to think more carefully about it. there is definitely something
     important about how it relates to the choice of dphi.

]]--
local calculate_beta = function(n_sub)
  local beta = {}

  local k = 1
  for i = 1, #n_sub do
    for j = 1, n_sub[i] do
      beta[k] = k
      k = k + 1
    end
  end
  beta[#beta+1] = k

  return beta
end

--[[ Interpolate theta using cubic spline interpolation.

    arguments:
    phi  - old grid coordinates.
    f    - value of function on old grid.
    nphi - new grid coordinates.

    returns:
    nf - f interpolated onto new grid coordinates.

    requires phi and nphi are monotonically increasing, and
    that #f == #phi.

--]]
local interp_cubic = function(phi, f, nphi)
  local b = 0
  local nf = {}

  assert(#phi == #f, "Function interp_cubic\n" ..
                     "Failed pre-condition #phi == #f\n" ..
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
    if not (b == #f) then
      local cdist = phi[b+1] - nphi[a]
      while cdist < 0 do
        b = b + 1
        if (b == #f) then break end
        cdist = phi[b+1] - nphi[a]
      end
    end

    local x, x0, x1, x2, x3
    local f0, f1, f2, f3
    if (b > 1) and (b < #f-1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b-1], phi[b], phi[b+1], phi[b+2]
      f0, f1, f2, f3 = f[b-1], f[b], f[b+1], f[b+2]
    elseif (b == 1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b], phi[b+1], phi[b+2], phi[b+3]
      f0, f1, f2, f3 = f[b], f[b+1], f[b+2], f[b+3]
    elseif (b == #f-1) then
      x, x0, x1, x2, x3 = nphi[a], phi[b-2], phi[b-1], phi[b], phi[b+1]
      f0, f1, f2, f3 = f[b-2], f[b-1], f[b], f[b+1]
    end

    if not ((b == #f) or (b == 0)) then
      nf[a] = f0*((x - x1)*(x - x2)*(x - x3))/((x0 - x1)*(x0 - x2)*(x0 - x3)) +
              f1*((x - x0)*(x - x2)*(x - x3))/((x1 - x0)*(x1 - x2)*(x1 - x3)) +
              f2*((x - x0)*(x - x1)*(x - x3))/((x2 - x0)*(x2 - x1)*(x2 - x3)) +
              f3*((x - x0)*(x - x1)*(x - x2))/((x3 - x0)*(x3 - x1)*(x3 - x2))
    elseif b == 0 then
      nf[a] = f[1]
    else
      nf[a] = f[#f]
    end

  end

  assert(#nf == #nphi, "Function interp_cubic\n" ..
                       "Failed post-condition #nphi == #nf\n" ..
                       "Exit\n\n") 

  return nf
end

--[[ Find peaks in free surface via linear interpolation

     finds zeros of theta by linear interpolation, and records location
     in terms of 'index' i into phi, and the value of phi(i)

     arguments:
     phi   - grid coordinates.
     v     - value of v on grid.
     u_mid - value of u on grid.

     returns:
     phi_peaks - value of phi at peaks
     i_peaks   - 'index' into grid for peak (non-integer)

     ensures that phi_peaks(i) lies in the correct interval of phi given i_peaks(i)
--]]
local find_peaks_uv = function(phi, v, u_mid)
  local phi_peaks = {}
  local i_peaks = {}

  assert(#phi == #v, "Function find_peaks\n" ..
                     "Failed pre-condition #phi == #theta\n" ..
                     "Exit\n\n")

  -- Experimental stuff
  -- need to know how 'peaky' a peak is
  -- Calculate the approximate measure of mean |curvature| for
  -- the entire surface.
  mean_curvature = 0
  for i = 2, (#v) do
    mean_curvature = mean_curvature + math.abs(v[i] - v[i-1]) /
                                      ((phi[i] - phi[i-1]) * math.abs(u_mid[i-1]))
  end
  mean_curvature = mean_curvature / (#v - 1)
  -- /Experimental stuff

  -- DODGY FIX i=2 -> i=3
  for i = 3, (#v) do
     if v[i] * v[i-1] < 0 and v[i-1] > 0 then
      -- We have found a 'peak' between i-1 and i, now we
      -- determine a more 'precise' location of the peak via linear interpolation.
      local t = v[i-1] / (v[i-1] - v[i])

      -- where did u go?
      -- Experimental stuff
      local dv = math.abs(v[i] - v[i-1])
      local dphi = phi[i] - phi[i-1]
      local um = math.abs(u_mid[i-1])
      if (dv / dphi) > mean_curvature then
      -- /Experimental stuff

        phi_peaks[#phi_peaks+1] = (1-t)*phi[i-1] + t*phi[i]
        i_peaks[#i_peaks+1]     = (i-1)+t

      else
	io.write('--- EXPERIMENTAL STUFF --- We are throwing out a peak!\n')
        print(i)
        io.flush()
      -- Experimental stuff
      end
      -- /Experimental stuff

    end
  end

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
--]]
local find_peaks = function(phi, theta, tau_mid)
  local phi_peaks = {}
  local i_peaks = {}

  assert(#phi == #theta, "Function find_peaks\n" ..
                         "Failed pre-condition #phi == #theta\n" ..
                         "Exit\n\n")

  mean_curvature = 0
  for i = 2, (#theta) do
    mean_curvature = mean_curvature + math.abs(theta[i] - theta[i-1]) / (phi[i] - phi[i-1])
  end
  mean_curvature = mean_curvature / (#theta - 1)

  -- DODGY FIX i=2 -> i=3
  for i = 3, (#theta) do
     if theta[i] * theta[i-1] < 0 and theta[i-1] > 0 then
      -- We have found a 'peak' between i-1 and i, now we
      -- determine a more 'precise' location of the peak via linear interpolation.
      local t = theta[i-1] / (theta[i-1] - theta[i])

      local dtheta = math.abs(theta[i] - theta[i-1])
      local dphi = phi[i] - phi[i-1]
      if (math.exp(tau_mid[i]) * dtheta / dphi) > mean_curvature then
 
        phi_peaks[#phi_peaks+1] = (1-t)*phi[i-1] + t*phi[i]
        i_peaks[#i_peaks+1]     = (i-1)+t

      else
	io.write('--- EXPERIMENTAL STUFF --- We are throwing out a peak!\n')
        io.flush()
  
      end

    end
  end
 
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

grid_functions.convert_grid_uv_hs = function(phi, v, u_mid, n, cluster_coeff, homotopy_s, phi_func)
  local compute_bsub = function(phi_crit_a, phi_crit_b, u_a, u_b)
     return (phi_crit_b - phi_crit_a) / math.pow(u_a + u_b,1)
  end

  local calculate_u = function(u_mid, k)
    local u_k = 0
    if k == 1 then
      u_k = (3*u_mid[k] - u_mid[k+1]) / 2
    elseif k == (#u_mid)+1 then
      u_k = (3*u_mid[k-1] - u_mid[k-2]) / 2
    else
      u_k = (u_mid[k-1] + u_mid[k]) / 2
    end
    
    return u_k
  end

  local phi_peaks, i_peaks = find_peaks_uv(phi, v, u_mid)

  local phi_crit = { 0 }
  for k, phi_peaks_k in pairs(phi_peaks) do
    phi_crit[k+1] = phi_peaks_k
  end
  phi_crit[#phi_crit+1] = phi[#phi]
  
  local i_crit = { 1 }
  for k, i_peaks_k in pairs(i_peaks) do
    i_crit[k+1] = i_peaks_k
  end
  i_crit[#i_crit+1] = #v

  local weighted_bsub_total = 0
  for k = 1, #phi_crit-1 do
    local uka = calculate_u(u_mid, math.floor(i_crit[k]))
    local ukb = calculate_u(u_mid, math.floor(i_crit[k+1]))

    uka = math.max(math.min(uka,1.0),0.0)
    ukb = math.max(math.min(ukb,1.0),0.0)

    local weighted_bsub_k = compute_bsub(phi_crit[k], phi_crit[k+1], uka, ukb)
    weighted_bsub_total = weighted_bsub_total + weighted_bsub_k
  end

  --[[
  Initial entries are related to the centre of the wave phi = 0
  --]]

  local u0a = calculate_u(u_mid, i_crit[1])
  u0a = math.max(math.min(u0a,1.0),0.0)

  local nn_sub = { }
  local nbeta_sub = { 1 }
  local nphi_sub = { phi_crit[1] }
  local ndphi_sub = { (cluster_coeff * math.pow(1-u0a,2) + (phi_crit[#phi_crit] * math.pow(u0a,2)))/n }

  --[[
    Next entries are determined by values at peaks
  --]]
  local cu_bsub_k = 0
  local cu_nsub_k = 0

  for k = 2, #phi_crit-1 do
    nphi_sub[k] = phi_crit[k]

    local uja = calculate_u(u_mid, math.floor(i_crit[k-1]))
    local ujb = calculate_u(u_mid, math.floor(i_crit[k]))
    uja = math.max(math.min(uja,1.0),0.0)
    ujb = math.max(math.min(ujb,1.0),0.0)
    ndphi_sub[k] = (cluster_coeff * math.pow(1-ujb,2) + (phi_crit[#phi_crit] * math.pow(ujb,2))) / n

    local del_bsub_k = compute_bsub(phi_crit[k-1], phi_crit[k], uja, ujb)
    
    cu_bsub_k = cu_bsub_k + del_bsub_k
    nn_sub[k-1] = math.floor(0.5 + (cu_bsub_k * (n-1) / weighted_bsub_total)) - cu_nsub_k
    cu_nsub_k = cu_nsub_k + nn_sub[k-1]
    nbeta_sub[k] = nbeta_sub[k-1] + nn_sub[k-1]
  end

  nphi_sub[#nphi_sub+1]   = phi_crit[#phi_crit]

  local uNa = calculate_u(u_mid, math.floor(i_crit[#i_crit-1]))
  local uNb = calculate_u(u_mid, math.floor(i_crit[#i_crit]))
  uNa = math.max(math.min(uNa,1.0),0.0)
  uNb = math.max(math.min(uNb,1.0),0.0)

  ndphi_sub[#ndphi_sub+1] = (cluster_coeff * math.pow(1-uNb,2) +  (phi_crit[#phi_crit] * math.pow(uNb,2))) / n
  local del_bsub_k = compute_bsub(phi_crit[#phi_crit-1], phi_crit[#phi_crit], uNa, uNb)

  cu_bsub_k = cu_bsub_k + del_bsub_k
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

  local beta = calculate_beta(nn_sub)
  local nphi = phi_func(beta, nbeta_sub, nphi_sub, ndphi_sub, homotopy_s)
  local nv   = interp_cubic(phi, v, nphi)

  return nv, nbeta_sub, nn_sub, nphi_sub, ndphi_sub

end

--[[ Converts grid from having no peak-data to having peak data. --]]
grid_functions.convert_grid_hs = function(phi, theta, tau_mid, n, cluster_coeff, homotopy_s, phi_func)

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

  local phi_peaks, i_peaks = find_peaks(phi, theta, tau_mid)

  --local phi_crit = { 0 }
  local phi_crit = { phi[1] }

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

    uka = math.max(math.min(uka,1.0),0.0)
    ukb = math.max(math.min(ukb,1.0),0.0)

    local weighted_bsub_k = compute_bsub(phi_crit[k], phi_crit[k+1], uka, ukb)
    weighted_bsub_total = weighted_bsub_total + weighted_bsub_k
  end

  --[[
  Initial entries are related to the centre of the wave phi = 0
  --]]

  local u0a = calculate_u(theta, tau_mid, i_crit[1])
  u0a = math.max(math.min(u0a,1.0),0.0)

  local nn_sub = { }
  local nbeta_sub = { 1 }
  local nphi_sub = { phi_crit[1] }
  local ndphi_sub = { (cluster_coeff * math.pow(1-u0a,2) + ((phi_crit[#phi_crit] - phi_crit[1]) * math.pow(u0a,2)))/n }

  --[[
    Next entries are determined by values at peaks
  --]]
  local cu_bsub_k = 0
  local cu_nsub_k = 0
  for k = 2, #phi_crit-1 do
    nphi_sub[k] = phi_crit[k]

    local uja = calculate_u(theta, tau_mid, math.floor(i_crit[k-1]))
    local ujb = calculate_u(theta, tau_mid, math.floor(i_crit[k]))
    uja = math.max(math.min(uja,1.0),0.0)
    ujb = math.max(math.min(ujb,1.0),0.0)
    ndphi_sub[k] = (cluster_coeff * math.pow(1-ujb,2) + ((phi_crit[#phi_crit] - phi_crit[1]) * math.pow(ujb,2))) / n

    local del_bsub_k = compute_bsub(phi_crit[k-1], phi_crit[k], uja, ujb)
    
    cu_bsub_k = cu_bsub_k + del_bsub_k
    nn_sub[k-1] = math.floor(0.5 + (cu_bsub_k * (n-1) / weighted_bsub_total)) - cu_nsub_k
    cu_nsub_k = cu_nsub_k + nn_sub[k-1]
    nbeta_sub[k] = nbeta_sub[k-1] + nn_sub[k-1]
  end

  nphi_sub[#nphi_sub+1]   = phi_crit[#phi_crit]

  local uNa = calculate_u(theta, tau_mid, math.floor(i_crit[#i_crit-1]))
  local uNb = calculate_u(theta, tau_mid, math.floor(i_crit[#i_crit]))
  uNa = math.max(math.min(uNa,1.0),0.0)
  uNb = math.max(math.min(uNb,1.0),0.0)

  ndphi_sub[#ndphi_sub+1] = (cluster_coeff * math.pow(1-uNb,2) +  ((phi_crit[#phi_crit] - phi_crit[1]) * math.pow(uNb,2))) / n
  local del_bsub_k = compute_bsub(phi_crit[#phi_crit-1], phi_crit[#phi_crit], uNa, uNb)

  cu_bsub_k = cu_bsub_k + del_bsub_k
  if #phi_crit > 2 then
    nn_sub[#nn_sub+1] = math.floor(0.5 + (cu_bsub_k * (n-1) / weighted_bsub_total)) - cu_nsub_k
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

  local beta   = calculate_beta(nn_sub)
  local nphi   = phi_func(beta, nbeta_sub, nphi_sub, ndphi_sub, homotopy_s)
  local ntheta = interp_cubic(phi, theta, nphi)

  return ntheta, nbeta_sub, nn_sub, nphi_sub, ndphi_sub

end

grid_functions.interp_cubic = interp_cubic
grid_functions.calculate_beta = calculate_beta

return grid_functions
