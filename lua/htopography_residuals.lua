require "surface"
require "integrals_sbox_homotopy"
require "grid_functions"

module("htopography_residuals", package.seeall)

bernoulli = function(tau, eta, froude)
  local f_sq  = froude * froude
  local alpha = math.exp(2.0*tau[#tau]) + (2.0*eta[#eta] / f_sq)
  local val = {}

  for i = 1, (#tau)-1 do
    val[i] = math.exp(2.0*tau[i]) +
             2.0 * (eta[i] / f_sq) -
             (alpha)
  end

  return val
end

compute_bttxe_hs = function(unpacked, x_calculation)

  -- Theta on free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of grid
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local n_sub      = unpacked.N_SUB
  local beta_sub   = unpacked.BETA_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  
  -- Parameters of topography
  local A     = unpacked.A
  local B     = unpacked.B
  local L     = unpacked.L
  local phi_c = unpacked.PHI_C

  -- Parameters of farstream
  local lambda     = unpacked.LAMBDA
  local D          = unpacked.D
  local gamma      = unpacked.GAMMA

  local n = #theta_s

  local beta_s = grid_functions.calculate_beta(n_sub)

  local tau_s_mid = integrals_sbox.tau_smoothbox_symmetric_s(theta_s, phi_sub,
                                                            dphi_sub,  beta_s,
                                                            beta_sub,       A,
                                                                   B,       L,
                                                               phi_c,       D,
                                                              lambda,   gamma,
                                                          homotopy_s)

  local theta_s_mid = {}
  for i = 1, n-1 do
     theta_s_mid[i] = (theta_s[i] + theta_s[i+1]) / 2.0
  end

  -- technically don't need this, except for output! should maybe include an option?
  local x_s_mid = {}
  local x_lhs = 0

  if x_calculation then
    x_s_mid, 
      x_lhs = integrals_sbox.x_symmetric_s(theta_s_mid, tau_s_mid,
                                                beta_s,  beta_sub,
                                               phi_sub,  dphi_sub, 
                                            homotopy_s)
  end

  local eta_s_mid,
          eta_lhs = integrals_sbox.y_symmetric_s(theta_s_mid, tau_s_mid,
                                                      beta_s,  beta_sub,
                                                     phi_sub,  dphi_sub,
                                                  homotopy_s)

  for j = 1, (n-1) do
    eta_s_mid[j] = (eta_s_mid[j] - eta_s_mid[#eta_s_mid])
  end

  return beta_s, theta_s_mid, tau_s_mid, x_s_mid, x_lhs, eta_s_mid, eta_lhs

end

--[[
compute_bttxe_hg = function(unpacked)

  -- Theta on free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of grid
  local n_sub      = unpacked.N_SUB
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S

  -- Parameters of topography
  local A = unpacked.A
  local B = unpacked.B
  local L = unpacked.L
  local phi_c = unpacked.PHI_C

  -- Parameters of far stream
  local lambda    = unpacked.LAMBDA
  local D_dstream = unpacked.D_DSTREAM
  local D_ustream = unpacked.D_USTREAM
  local gamma     = unpacked.GAMMA

  local n = #theta_s

  local beta_s = grid_functions.calculate_beta(n_sub)

  -- need homotopy parameter etc.
  local tau_s_mid = integrals_sbox.tau_smoothbox_s(theta_s,    phi_sub,
                                                  dphi_sub,     beta_s,
                                                  beta_sub,          A,
                                                         B,          L,
                                                     phi_c,  D_dstream,
                                                 D_ustream,     lambda,
                                                     gamma, homotopy_s)

  local theta_s_mid = {}
  for i = 1, n-1 do
     theta_s_mid[i] = (theta_s[i] + theta_s[i+1]) / 2.0
  end

  -- technically don't need this, except for output! should maybe include an option?
  local x_s_mid = {}
  local x_lhs = 0

  if x_calculation then
    x_s_mid, 
      x_lhs = integrals_sbox.x_s(theta_s_mid, tau_s_mid,
                                      beta_s,  beta_sub,
                                     phi_sub,  dphi_sub, 
                                   homotopy_s)
  end

  local eta_s_mid,
          eta_lhs = integrals_sbox.y_s(theta_s_mid, tau_s_mid,
                                            beta_s,  beta_sub,
                                           phi_sub,  dphi_sub,
                                        homotopy_s)

  for j = 1, (n-1) do
    eta_s_mid[j] = (eta_s_mid[j] - eta_s_mid[#eta_s_mid])
  end

  return beta_s, theta_s_mid, tau_s_mid, x_s_mid, x_lhs, eta_s_mid, eta_lhs
  
end
--]]

--[[
compute_output_hg = function(unpacked)

  local t = {}
    
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hg(unpacked, true)

  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end

  for i = 1, #eta_s_mid do
    eta_s_mid[i] = (eta_s_mid[i] - eta_s_mid[#eta_s_mid])
  end

  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  local phi_s_mid  = integrals_sbox.phi(beta_s_mid,
                                          beta_sub,
                                           phi_sub,
                                          dphi_sub,
                                        homotopy_s)

  -- Free surface variables
  t.theta_s_mid = theta_s_mid
  t.tau_s_mid   = tau_s_mid
  t.x_s_mid     = x_s_mid
  t.eta_s_mid   = eta_s_mid
  t.phi_s_mid   = phi_s_mid

  -- Parameters of flow/topography
  t.A     = unpacked.A
  t.B     = unpacked.B
  t.L     = unpacked.L
  t.phi_c = unpacked.PHI_C
  t.F     = unpacked.FROUDE

  -- Parameters of far stream
  t.D_dstream = unpacked.D_DSTREAM
  t.D_ustream = unpacked.D_USTREAM
  t.lambda    = unpacked.LAMBDA
  t.gamma     = unpacked.GAMMA

  -- Parameters of grid
  t.n_sub      = unpacked.N_SUB
  t.phi_sub    = unpacked.PHI_SUB
  t.dphi_sub   = unpacked.DPHI_SUB
  t.homotopy_s = unpacked.HOMOTOPY_S
  
  return t
  
end
--]]

compute_output_hs = function(unpacked)

  local t = {}
    
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hs(unpacked, true)

  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end

  for i = 1, #eta_s_mid do
    eta_s_mid[i] = (eta_s_mid[i] - eta_s_mid[#eta_s_mid])
  end

  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  local phi_s_mid  = integrals_sbox.phi(beta_s_mid,
                                          beta_sub,
                                           phi_sub,
                                          dphi_sub,
                                        homotopy_s)

  -- Hmm!
  local theta_s = unpacked.THETA_S
  local A       = unpacked.A
  local B       = unpacked.B
  local L       = unpacked.L
  local D       = unpacked.D
  local phi_c   = unpacked.PHI_C
  local lambda  = unpacked.LAMBDA
  local gamma   = unpacked.GAMMA
  
  local topography_start = os.clock()
  
  local x_b_mid = integrals_sbox.x_smoothbox_symmetric_b(phi_s_mid,    theta_s,
                                                            beta_s,    phi_sub,
                                                          dphi_sub,   beta_sub,
                                                                 A,          B,
                                                                 L,      phi_c,
                                                                 D,     lambda,
								 gamma, homotopy_s)
  local y_b_mid = integrals_sbox.y_smoothbox_symmetric_b(phi_s_mid,    theta_s,
                                                            beta_s,    phi_sub,
                                                          dphi_sub,   beta_sub,
                                                                 A,          B,
                                                                 L,      phi_c,
                                                                 D,     lambda,
                                                             gamma, homotopy_s)

  for i = 1, #y_b_mid do
    y_b_mid[i] = y_b_mid[i] - y_b_mid[#y_b_mid]
  end
  local topography_finish = os.clock()
  io.write('    Topography time = ' .. tostring(topography_finish - topography_start) ..'\n')
  io.flush()
  -- Free surface variables
  t.theta_s_mid = theta_s_mid
  t.tau_s_mid   = tau_s_mid
  t.x_s_mid     = x_s_mid
  t.eta_s_mid   = eta_s_mid
  t.phi_s_mid   = phi_s_mid

  -- hmm
  t.x_b_mid = x_b_mid
  t.y_b_mid = y_b_mid

  -- Parameters of flow/topography
  t.A      = A
  t.B      = B
  t.L      = L
  t.phi_c  = phi_c
  t.FROUDE = unpacked.FROUDE
  if not (unpacked.ETA0 == nil) then
    t.eta0   = unpacked.ETA0
  end
  if not (unpacked.U_MIN == nil) then
    t.u_min = unpacked.U_MIN
  end

  -- Parameters of far stream
  t.D      = D
  t.lambda = lambda
  t.gamma  = gamma

  -- Parameters of grid
  t.n_sub       = unpacked.N_SUB
  t.beta_sub    = unpacked.BETA_SUB
  t.phi_sub     = phi_sub
  t.dphi_sub    = dphi_sub
  t.homotopy_s  = homotopy_s
  t.cluster_val = unpacked.CLUSTER_VAL
  
  return t
  
end

unit_clamp = function(g)
  return math.max(math.min(g,1), 0)
end

residual_hs = function(unknowns, extra)
  local fs = extra[1]

  if print_out then
    io.write('htopography_residuals.residual_hs entry point ... \n')
    io.flush()
  end

  local unpacked = surface.unpack_vec(unknowns,
                                    fs.initial,
                                    fs.is_free,
                               fs.is_continued,
                                  fs.ind_table)

  -- Theta on the free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of topography
  local A     = unpacked.A
  local B     = unpacked.B
  local L     = unpacked.L
  local phi_c = unpacked.PHI_C

  -- Flow parameters
  local froude     = unpacked.FROUDE
  local eta0       = unpacked.ETA0

  -- Far-stream parameters
  local lambda = unpacked.LAMBDA
  local D      = unpacked.D
  local gamma  = unpacked.GAMMA

  -- Grid parameters
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local n_sub      = unpacked.N_SUB
  local n_b        = unpacked.N_B
  local homotopy_s = unpacked.HOMOTOPY_S

  -- Optional flow parameters
  local u_min       = unpacked.U_MIN
  local cluster_val = unpacked.CLUSTER_VAL

  if print_out then
    io.write('About to compute beta, tau, theta etc ...')
  end

  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hs(unpacked)
  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end
  if print_out then
    io.write(' done\n.')
    io.flush()
  end
 
  local phi_s_mid  = integrals_sbox.phi(beta_s_mid,
                                          beta_sub,
                                           phi_sub,
                                          dphi_sub,
                                        homotopy_s)

  if print_out then
    io.flush()
    io.write('n_sub = [')
    for kk,vv in pairs(n_sub) do
      io.write(tostring(vv) .. ' ')
    end
    io.write('];\n')
    io.flush()

    io.write('\nphi_sub = [')
    for kk,vv in pairs(phi_sub) do
      io.write(tostring(vv) .. ' ')
    end
    io.write('];\n')
    io.flush()
    
    io.write('dphi_sub = [')
    for kk,vv in pairs(dphi_sub) do
      io.write(tostring(vv) .. ' ')
    end
    io.write('];\n')
    io.flush()
  end

  local n = (#theta_s)
  if print_out then
    io.write('About to compute bernoullis eqn ...')
  end

  -- N-2 equations
  local val = bernoulli(tau_s_mid, eta_s_mid, froude)
  if print_out then
    io.write(' done.\n')
    io.write('Compute crest equations and linearised stuff ...')
    io.flush()
  end 

  -- (N-2) + #crests
  -- Equations for crest(s)
  local k = 0
  for i = 1, (#n_sub-1) do
     k = k + n_sub[i]
     val[#val+1] = (theta_s_mid[k] + theta_s_mid[k+1])
     --val[#val+1] = beta_sub[i+1] - (k+1)
  end

  -- NEED TO DOUBLE CHECK FF EQUATIONS (looks like they make sense if D<0)
  -- (N) + #crests
  -- Far field matching equations
  local phi_s_ff_a   = phi_s_mid[n-1]
  local u_s_ff_a     = 1 + (D*math.cos(lambda)*math.exp(-lambda*phi_s_ff_a))
  local v_s_ff_a     =     (D*math.sin(lambda)*math.exp(-lambda*phi_s_ff_a))
  local theta_s_ff_a = math.atan(v_s_ff_a / u_s_ff_a)
  --local tau_s_ff_a   = math.log(u_s_ff_a * u_s_ff_a + v_s_ff_a * v_s_ff_a) / 2.0

  val[#val+1] = theta_s_mid[n-1] - theta_s_ff_a
  --val[#val+1] = tau_s_mid[n-1] - tau_s_ff_a

  phi_s_ff_b  = phi_s_mid[n-2]
  u_s_ff_b    = 1 + (D*math.cos(lambda)*math.exp(-lambda*phi_s_ff_b))
  v_s_ff_b    =     (D*math.sin(lambda)*math.exp(-lambda*phi_s_ff_b))
  theta_s_ff_b = math.atan(v_s_ff_b / u_s_ff_b)
  val[#val+1] = theta_s_mid[n-2] - theta_s_ff_b

  -- (N+1) + #crests
  -- Equation for linearised far field
  val[#val+1] = math.tan(lambda) - froude*froude*lambda

  -- (N+2) + #crests
  -- Geometry equations
  local x_b = integrals_sbox.x_smoothbox_symmetric_b({phi_c},    theta_s,
                                                      beta_s,    phi_sub,
                                                    dphi_sub,   beta_sub,
                                                           A,          B,
                                                           L,      phi_c,
                                                           D,     lambda,
                                                       gamma, homotopy_s)
  val[#val+1] = x_b[1] - (L/2.0)
  val[#val+1] = eta_s_mid[1] - eta0

  if print_out then
    io.write(' done.\n')
    io.write('Optional parameters ... ')
    io.flush()
  end

  -- Optional parameters
  -- (N+2) + #crests + u_min
 if u_min then
    local integral_u_min = 1/0
    for kk, vv in pairs(tau_s_mid) do
      if (math.cos(theta_s_mid[kk])*math.exp(vv)) < integral_u_min then
        integral_u_min = math.cos(theta_s_mid[kk])*math.exp(vv)
      end
    end
    val[#val+1] = u_min - integral_u_min
  end

  -- (N+3) + #crests + u_min + #n_sub
  if cluster_val then
    local k = 1
    for i = 1, (#n_sub) do
      local u_k = unit_clamp(math.cos(theta_s_mid[k])*math.exp(tau_s_mid[k]))
      if print_out_box then print(u_k) end
      val[#val+1] = dphi_sub[i] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n)
      k = k + n_sub[i]
    end
    assert(k == n, 'Failed assertion k == n for some thing')
    local u_k = unit_clamp(math.cos(theta_s_mid[#theta_s_mid])*math.exp(tau_s_mid[#tau_s_mid]))
    if print_out_box then print(u_k) end
    val[#val+1] = dphi_sub[#dphi_sub] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n)
  end
  
  if print_out then
    io.write(' done.\nValue of residuals :\n')
    for kk, vv in pairs(val) do
      io.write(tostring(kk) .. ' : ' .. tostring(vv) .. '.\n')
    end
    io.flush()
  end
  
  return val
end

-- **************************************************************************************************************** --
-- **************************************************************************************************************** --
--[[
residual_hg = function(unknowns, extra)
  local fs = extra[1]

  local unpacked = surface.unpack_vec(unknowns,
                                    fs.initial,
                                    fs.is_free,
                               fs.is_continued,
                                  fs.ind_table)
                                  
  -- Theta on the free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of topography
  local topography_a = unpacked.TOPOGRAPHY_A
  local topography_b = unpacked.TOPOGRAPHY_B

  -- Flow parameters
  local froude     = unpacked.FROUDE
  local eta0       = unpacked.ETA0

  -- Far-stream parameters
  local lambda     = unpacked.LAMBDA
  local Du_dstream = unpacked.Du_DSTREAM
  local Dv_dstream = unpacked.Dv_DSTREAM
  local Du_ustream = unpacked.Du_USTREAM
  local Dv_ustream = unpacked.Dv_USTREAM
  local gamma      = unpacked.GAMMA

  -- Grid parameters
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local n_sub      = unpacked.N_SUB
  local n_b        = unpacked.N_B
  local homotopy_s = unpacked.HOMOTOPY_S

  -- Optional flow parameters
  local u_min       = unpacked.U_MIN
  local cluster_val = unpacked.CLUSTER_VAL

  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hg(unpacked)
  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end
  
  local phi_s_mid  = integrals_sbox.phi(beta_s_mid,
                                          beta_sub,
                                           phi_sub,
                                          dphi_sub,
                                        homotopy_s)

  if print_out then
    io.flush()
    io.write('\nphi_sub = [')
    for kk,vv in pairs(phi_sub) do
      io.write(tostring(vv) .. ' ')
    end
    io.write('];\n')
    io.flush()
    
    io.write('dphi_sub = [')
    for kk,vv in pairs(dphi_sub) do
      io.write(tostring(vv) .. ' ')
    end
    io.write('];\n')
    io.flush()
  end

  local n = (#theta_s)

  if print_out then
    io.write('About to compute bernoullis eqn ...')
  end

  local val = bernoulli(tau_s_mid, eta_s_mid, froude) -- N-3 equations
  if print_out then
    io.write(' done.\n')
    io.write('Compute crest equations and linearised stuff ...')
    io.flush()
  end 

  -- Equations for crest(s)
  local k = 0
  for i = 1, (#n_sub-1) do
     k = k + n_sub[i]
     val[#val+1] = (theta_s_mid[k] + theta_s_mid[k+1])
  end

  -- Far field matching equations
  local phi_asymp_a = phi_s_mid[n-1]
  local u_lim_asymp_a   = 1 - (Du_dstream*math.cos(lambda)*math.exp(-lambda*phi_asymp_a))
  local v_lim_asymp_a   = -(Dv_dstream*math.sin(lambda)*math.exp(-lambda*phi_asymp_a))

  local tau_lim_asymp_a = math.log(u_lim_asymp_a * u_lim_asymp_a + v_lim_asymp_a * v_lim_asymp_a) / 2.0
  val[#val+1] = tau_s_mid[n-1] - tau_lim_asymp_a

  local theta_lim_asymp_a = math.atan(v_lim_asymp_a / u_lim_asymp_a)
  val[#val+1] = theta_s_mid[n-1] - theta_lim_asymp_a

  local phi_asymp_c = phi_s_mid[1]
  local u_lim_asymp_c   = 1 - (Du_ustream*math.cos(lambda)*math.exp(lambda*phi_asymp_c))
  local v_lim_asymp_c   = -(Dv_ustream*math.sin(lambda)*math.exp(lambda*phi_asymp_c))

  local tau_lim_asymp_c = math.log(u_lim_asymp_c * u_lim_asymp_c + v_lim_asymp_c * v_lim_asymp_c) / 2.0
  val[#val+1] = tau_s_mid[1] - tau_lim_asymp_c

  -- Dv_dstream =  Du_dstream,  Dv_ustream = -Du_ustream
  val[#val+1] = Dv_dstream - Du_dstream
  val[#val+1] = Dv_ustream + Du_ustream

  local theta_lim_asymp_c = math.atan(v_lim_asymp_c / u_lim_asymp_c)
  val[#val+1] = theta_s_mid[1] - theta_lim_asymp_c

  -- Equation for linearised far field
  val[#val+1] = math.tan(lambda) - froude*froude*lambda -- 122

  -- Geometry equations
  val[#val+1] = eta_s_mid[math.floor(n/2)] - eta0 -- 127

  if print_out then
    io.write(' done.\n')
    io.write('Optional parameters ... ')
    io.flush()
  end
    
  -- Optional parameters
  if u_min then
    local integral_u_min = 1/0
    for kk, vv in pairs(tau_s_mid) do
      if (math.cos(theta_s_mid[kk])*math.exp(vv)) < integral_u_min then
        integral_u_min = math.cos(theta_s_mid[kk])*math.exp(vv)
      end
    end
    val[#val+1] = u_min - integral_u_min
  end

  if cluster_val then
    local k = 1
    for i = 1, (#n_sub) do
      local u_k = unit_clamp(math.cos(theta_s_mid[k])*math.exp(tau_s_mid[k]))
      if print_out_box then print(u_k) end
      val[#val+1] = dphi_sub[i] - (((cluster_val * math.pow(1-u_k,2)) +
                                   ((phi_sub[#phi_sub] - phi_sub[1])  * math.pow(u_k,2)))/n)
      k = k + n_sub[i]
    end
    assert(k == n, 'Failed assertion k == n for some thing')
    local u_k = unit_clamp(math.cos(theta_s_mid[#theta_s_mid])*math.exp(tau_s_mid[#tau_s_mid]))
    if print_out_box then print(u_k) end
    val[#val+1] = dphi_sub[#dphi_sub] - (((cluster_val * math.pow(1-u_k,2)) +
                                         ((phi_sub[#phi_sub] - phi_sub[1]) * math.pow(u_k,2)))/n)
  end
  
  if print_out then
    io.write(' done.\nValue of residuals :\n')
    for kk, vv in pairs(val) do
      io.write(tostring(kk) .. ' : ' .. tostring(vv) .. '.\n')
    end
    io.flush()
  end

  return val

end
--]]

