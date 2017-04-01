require "surface"
require "integrals_box_homotopy"

module("sbox_residuals", package.seeall)


--[[ Calculate value of beta for grid

     beta has been reduced to a cumulative count of how many grid points
     there are. this seems redundant as it must act like an index, but
     i'd have to think more carefully about it. there is definitely something
     important about how it relates to the choice of dphi.

]]--
calculate_beta = function(n_sub)
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

bernoulli = function(tau, eta, froude)
  local f_sq  = froude * froude
  local alpha = math.exp(2*tau[#tau])
  local val = {}

  for i = 1, (#tau)-1 do
    val[i] =  math.exp(2*tau[i]) +
              2.0 * ((eta[i] - eta[#eta]) / f_sq) -
              (alpha)
  end

  return val
end

compute_bttxe_hs = function(unpacked)

  -- Theta on free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of farstream
  local lambda     = unpacked.LAMBDA
  local Du         = unpacked.Du
  local Dv         = unpacked.Dv
  local gamma      = unpacked.GAMMA

  -- Parameters of box
  local phi_box    = unpacked.PHI_BOX
  local box_height = unpacked.BOX_HEIGHT

  -- Parameters of grid
  local n_sub      = unpacked.N_SUB
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S

  local n = #theta_s

  local beta_s = calculate_beta(n_sub)

  local tau_s_mid = integrals_homotopy.tau_symmetric_s(theta_s, phi_sub, dphi_sub,
                                                       beta_s, beta_sub, phi_box, Du, Dv,
                                                       lambda, box_height, gamma,
                                                       homotopy_s)

  local theta_s_mid = {}
  for i = 1, n-1 do
     theta_s_mid[i] = (theta_s[i] + theta_s[i+1]) / 2.0
  end

  local x_s_mid, 
          x_lhs = integrals_homotopy.x_symmetric_s(theta_s_mid, tau_s_mid,
                                                  beta_s, beta_sub, phi_sub, dphi_sub, 
                                                  homotopy_s)
  local eta_s_mid,
          eta_lhs = integrals_homotopy.y_symmetric_s(theta_s_mid, tau_s_mid,
                                                     beta_s, beta_sub, phi_sub, dphi_sub,
                                                     homotopy_s)

  for j = 1, (n-1) do
    eta_s_mid[j] = (eta_s_mid[j] - eta_s_mid[#eta_s_mid])
  end

  return beta_s, theta_s_mid, tau_s_mid, x_s_mid, x_lhs, eta_s_mid, eta_lhs

end

compute_bttxe_hg = function(unpacked)

  -- Theta on free surface
  local theta_s    = unpacked.THETA_S

  -- Parameters of far stream
  local lambda     = unpacked.LAMBDA
  local Du_dstream = unpacked.Du_DSTREAM
  local Dv_dstream = unpacked.Dv_DSTREAM
  local Du_ustream = unpacked.Du_USTREAM
  local Dv_ustream = unpacked.Dv_USTREAM
  local gamma      = unpacked.GAMMA

  -- Parameters of box
  local phi_box    = unpacked.PHI_BOX
  local box_height = unpacked.BOX_HEIGHT

  -- Parameters of grid
  local n_sub      = unpacked.N_SUB
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S

  local n = #theta_s

  local beta_s = calculate_beta(n_sub)

  -- need homotopy parameter etc.
  local tau_s_mid = integrals_homotopy.tau_s(theta_s, phi_sub, dphi_sub,
                                             beta_s, beta_sub, phi_box, 
                                             Du_dstream, Dv_dstream,
                                             Du_ustream, Dv_ustream,
                                             lambda, box_height, gamma,
                                             homotopy_s)

  local theta_s_mid = {}
  for i = 1, n-1 do
     theta_s_mid[i] = (theta_s[i] + theta_s[i+1]) / 2.0
  end

  local x_s_mid, 
          x_lhs = integrals_homotopy.x_s(theta_s_mid, tau_s_mid,
                                         beta_s, beta_sub, phi_sub, dphi_sub, 
                                         homotopy_s)
  local eta_s_mid,
          eta_lhs = integrals_homotopy.y_s(theta_s_mid, tau_s_mid,
                                           beta_s, beta_sub, phi_sub, dphi_sub,
                                           homotopy_s)
  
  for j = 1, (n-1) do
    eta_s_mid[j] = (eta_s_mid[j] - eta_s_mid[#eta_s_mid])
  end

  return beta_s, theta_s_mid, tau_s_mid, x_s_mid, x_lhs, eta_s_mid, eta_lhs
  
end

compute_output_hg = function(unpacked)

  local t = {}
    
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_homotopy(unpacked)

  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end

  for i = 1, #eta_s_mid do
    eta_s_mid[i] = (eta_s_mid[i] - eta_s_mid[#eta_s_mid])
  end

  local phi_box    = unpacked.PHI_BOX
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  local phi_s_mid  = integrals_homotopy.phi(beta_s_mid, beta_sub, phi_sub, dphi_sub, homotopy_s)

  -- Free surface variables
  t.theta_s_mid = theta_s_mid
  t.tau_s_mid   = tau_s_mid
  t.x_s_mid     = x_s_mid
  t.eta_s_mid   = eta_s_mid
  t.phi_s_mid   = phi_s_mid

  -- Parameters of flow
  t.phi_box    = unpacked.PHI_BOX
  t.box_height = unpacked.BOX_HEIGHT
  t.box_width  = unpacked.BOX_WIDTH
  t.F          = unpacked.FROUDE

  -- Parameters of far stream
  t.Du_dstream = unpacked.Du_DSTREAM
  t.Dv_dstream = unpacked.Dv_DSTREAM
  t.Du_ustream = unpacked.Du_USTREAM
  t.Dv_ustream = unpacked.Dv_USTREAM
  t.lambda     = unpacked.LAMBDA
  t.gamma      = unpacked.GAMMA

  -- Parameters of grid
  t.n_sub      = unpacked.N_SUB
  t.phi_sub    = unpacked.PHI_SUB
  t.dphi_sub   = unpacked.DPHI_SUB
  t.homotopy_s = unpacked.HOMOTOPY_S
  
  return t
  
end

compute_output_hs = function(unpacked)

  local t = {}
    
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hs(unpacked)

  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end

  for i = 1, #eta_s_mid do
    eta_s_mid[i] = (eta_s_mid[i] - eta_s_mid[#eta_s_mid])
  end

  local phi_box    = unpacked.PHI_BOX
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  local phi_s_mid  = integrals_homotopy.phi(beta_s_mid, beta_sub, phi_sub, dphi_sub, homotopy_s)

  -- Free surface variables
  t.theta_s_mid = theta_s_mid
  t.tau_s_mid   = tau_s_mid
  t.x_s_mid     = x_s_mid
  t.eta_s_mid   = eta_s_mid
  t.phi_s_mid   = phi_s_mid

  -- Parameters of flow
  t.phi_box    = unpacked.PHI_BOX
  t.box_height = unpacked.BOX_HEIGHT
  t.box_width  = unpacked.BOX_WIDTH
  t.F          = unpacked.FROUDE

  -- Parameters of far stream
  t.Du     = unpacked.Du
  t.Dv     = unpacked.Dv
  t.lambda = unpacked.LAMBDA
  t.gamma  = unpacked.GAMMA

  -- Parameters of grid
  t.n_sub      = unpacked.N_SUB
  t.phi_sub    = unpacked.PHI_SUB
  t.dphi_sub   = unpacked.DPHI_SUB
  t.homotopy_s = unpacked.HOMOTOPY_S
  
  return t
  
end

unit_clamp = function(g)
  return math.max(math.min(g,1), 0)
end

residual_hs = function(unknowns, extra)
  local fs = extra[1]
  if print_out_box then
    print('here_6_1')
  end
  local unpacked = surface.unpack_vec(unknowns,
                                    fs.initial,
                                    fs.is_free,
                               fs.is_continued,
                                  fs.ind_table)

  -- Theta on the free surface
  local theta_s    = unpacked.THETA_S

  -- Box parameters
  local phi_box    = unpacked.PHI_BOX
  local box_height = unpacked.BOX_HEIGHT
  local box_width  = unpacked.BOX_WIDTH

  -- Flow parameters
  local froude     = unpacked.FROUDE
  local eta0       = unpacked.ETA0

  -- Far-stream parameters
  local lambda     = unpacked.LAMBDA
  local Du         = unpacked.Du
  local Dv         = unpacked.Dv
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

  if print_out_box then
    print('here_6_2')
  end
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs,
                  eta_s_mid,   eta_lhs = compute_bttxe_hs(unpacked)
  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end
  if print_out_box then
    print('here_6_3')
  end
  local phi_s_mid  = integrals_homotopy.phi(beta_s_mid, beta_sub, phi_sub, dphi_sub, homotopy_s)
  if print_out_box then
  for kk,vv in pairs(phi_sub) do
     print(kk,vv)
  end
  for kk,vv in pairs(dphi_sub) do
     print(kk,vv)
  end
  end

  local n = (#theta_s)

  local phi_b_03 = {}
  local phi_b_34 = {}
    
  for i = 1, n_b do
    phi_b_03[i] = ((i-1)*phi_box[3])/(n_b - 1)
    phi_b_34[i] = ((n_b-i)*phi_box[3] + (i-1)*phi_box[4])/(n_b - 1)
  end

  local phi_b_03_mid = {}
  local phi_b_34_mid = {}
  for i = 1, n_b-1 do
    phi_b_03_mid[i] = (phi_b_03[i] + phi_b_03[i+1]) / 2.0
    phi_b_34_mid[i] = (phi_b_34[i] + phi_b_34[i+1]) / 2.0
  end

  if print_out_box then
    print('here_6_4')
    io.flush()
  end

  local tau_b_03_mid = integrals_homotopy.tau_symmetric_b(theta_s, phi_sub, dphi_sub,
                                                          beta_s, beta_sub, phi_box, phi_b_03_mid,
                                                          Du, Dv, lambda, box_height, gamma,
                                                          homotopy_s)
  if print_out_box then
    print('here_6_5')
  end
  local tau_b_34_mid = integrals_homotopy.tau_symmetric_b(theta_s, phi_sub, dphi_sub,
                                                          beta_s, beta_sub, phi_box, phi_b_34_mid,
                                                          Du, Dv, lambda, box_height, gamma,
                                                          homotopy_s)
  if print_out_box then
    print('here_6_6')
  end
  local x_03 = integrals_homotopy.geometry_b(tau_b_03_mid,phi_b_03)
  if print_out_box then
    print('here_6_7')
  end
  local y_34 = integrals_homotopy.geometry_b(tau_b_34_mid,phi_b_34)

  local val = bernoulli(tau_s_mid, eta_s_mid, froude) -- N-2 equations

  -- Equations for crest(s)
  local k = 0
  for i = 1, (#n_sub-1) do
     k = k + n_sub[i]
     val[#val+1] = (theta_s_mid[k] + theta_s_mid[k+1])
     --val[#val+1] = beta_sub[i+1] - (k+1)
  end

  -- Far field matching equations
  u_lim_asymp     = 1 - (Du*math.cos(lambda)*math.exp(-lambda*phi_s_mid[n-1]))
  v_lim_asymp     = -(Dv*math.sin(lambda)*math.exp(-lambda*phi_s_mid[n-1]))
  tau_lim_asymp   = math.log(u_lim_asymp * u_lim_asymp + v_lim_asymp * v_lim_asymp) / 2.0

  val[#val+1] = tau_s_mid[n-1] - tau_lim_asymp -- 119

  --phi_asymp_b = integrals_cluster.phi({beta[n]}, phi_sub, dphi_sub)
  --phi_asymp_b = phi_asymp_b[1]
  phi_asymp_b = phi_s_mid[n-1]
  u_lim_asymp     = 1 - (Du*math.cos(lambda)*math.exp(-lambda*phi_asymp_b))
  v_lim_asymp     = -(Dv*math.sin(lambda)*math.exp(-lambda*phi_asymp_b))
  theta_lim_asymp_b = math.atan(v_lim_asymp / u_lim_asymp)

  val[#val+1] = theta_s_mid[n-1] - theta_lim_asymp_b -- 120

  --phi_asymp_a = integrals_cluster.phi({beta[n-1]}, phi_sub, dphi_sub)
  --phi_asymp_a = phi_asymp_a[1]
  phi_asymp_a = phi_s_mid[n-2]

  u_lim_asymp     = 1 - (Du*math.cos(lambda)*math.exp(-lambda*phi_asymp_a))
  v_lim_asymp     = -(Dv*math.sin(lambda)*math.exp(-lambda*phi_asymp_a))
  theta_lim_asymp_a = math.atan(v_lim_asymp / u_lim_asymp)

  val[#val+1] = theta_s_mid[n-2] - theta_lim_asymp_a -- 121

  -- Equation for linearised far field
  val[#val+1] = math.tan(lambda) - froude*froude*lambda -- 122

  -- Geometry equations
  val[#val+1] = x_03 - (box_width/2.0) -- 123
  val[#val+1] = y_34 - math.abs(box_height) -- 124
  val[#val+1] = phi_box[2] + phi_box[3] -- 125
  val[#val+1] = phi_box[1] + phi_box[4] -- 126
  val[#val+1] = eta_s_mid[1] - eta0 -- 127

  if print_out_box then
    print('here_6_8')
  end
  -- Optional parameters
  if u_min then
    local integral_u_min = 1/0
    for kk, vv in pairs(tau_s_mid) do
      if (math.cos(theta_s_mid[kk])*math.exp(vv)) < integral_u_min then
        integral_u_min = math.cos(theta_s_mid[kk])*math.exp(vv)
      end
    end
    val[#val+1] = u_min - integral_u_min -- 128
  end

  if print_out_box then
    for kk, vv in pairs(phi_s_mid) do
      print(kk, vv)
    end
  end

  if print_out_box then
    print('here_6_9')
  end
  if cluster_val then
    local k = 1
    for i = 1, (#n_sub) do
      local u_k = unit_clamp(math.cos(theta_s_mid[k])*math.exp(tau_s_mid[k]))
      if print_out_box then print(u_k) end
      val[#val+1] = dphi_sub[i] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n) -- 129
      k = k + n_sub[i]
    end
    assert(k == n, 'Failed assertion k == n for some thing')
    local u_k = unit_clamp(math.cos(theta_s_mid[#theta_s_mid])*math.exp(tau_s_mid[#tau_s_mid]))
    if print_out_box then print(u_k) end
    val[#val+1] = dphi_sub[#dphi_sub] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n) -- 130
  end

  if print_out_box then
    for kk, vv in pairs(val) do
      print(kk, vv)
    end
    io.flush()
  end


  return val
end

residual_hg = function(unknowns, extra)
  local fs = extra[1]
  
  if print_out_box then
    print('here_6_1')
  end
  local unpacked = surface.unpack_vec(unknowns,
                                    fs.initial,
                                    fs.is_free,
                               fs.is_continued,
                                  fs.ind_table)
  local box_height = unpacked.BOX_HEIGHT
  local box_width  = unpacked.BOX_WIDTH
  local froude     = unpacked.FROUDE
  local lambda     = unpacked.LAMBDA
  local Du_dstream = unpacked.Du_DSTREAM
  local Dv_dstream = unpacked.Dv_DSTREAM
  local Du_ustream = unpacked.Du_USTREAM
  local Dv_ustream = unpacked.Dv_USTREAM
  local n_sub      = unpacked.N_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local theta_s    = unpacked.THETA_S
  local n_b        = unpacked.N_B
  local phi_box    = unpacked.PHI_BOX
  local gamma      = unpacked.GAMMA

  local f_sq = froude * froude

  if print_out_box then
    print('here_6_1')
  end
  local beta_s, theta_s_mid, tau_s_mid,
                    x_s_mid,     x_lhs, 
                  eta_s_mid,   eta_lhs = compute_ttxep(unpacked)

  local m  = #n_sub
  local n = (#theta_s)

  local phi_b_12 = {}
  local phi_b_20 = {}
  local phi_b_03 = {}
  local phi_b_34 = {}
  for i = 1, n_b do
    phi_b_12[i] = ((n_b-i)*phi_box[1] + (i-1)*phi_box[2])/(n_b - 1)
    phi_b_20[i] = ((n_b-i)*phi_box[2])/(n_b - 1)
    phi_b_03[i] = ((i-1)*phi_box[3])/(n_b - 1)
    phi_b_34[i] = ((n_b-i)*phi_box[3] + (i-1)*phi_box[4])/(n_b - 1)
  end

  local phi_b_12_mid = {}
  local phi_b_20_mid = {}
  local phi_b_03_mid = {}
  local phi_b_34_mid = {}
  for i = 1, n_b-1 do
    phi_b_12_mid[i] = (phi_b_12[i] + phi_b_12[i+1]) / 2.0
    phi_b_20_mid[i] = (phi_b_20[i] + phi_b_20[i+1]) / 2.0
    phi_b_03_mid[i] = (phi_b_03[i] + phi_b_03[i+1]) / 2.0
    phi_b_34_mid[i] = (phi_b_34[i] + phi_b_34[i+1]) / 2.0
  end

  if print_out_box then
    print('here_6_1')
  end


  local tau_b_12_mid = integrals.tau_cluster_b(theta_s, phi_sub, dphi_sub,
                                               beta_s, beta_sub, phi_box, phi_b_12_mid,
                                               Du_dstream, Dv_dstream, 
                                               Du_ustream, Dv_ustream,
                                               lambda, box_height, gamma)
  if print_out_box then
    print('here_6_1')
  end
  local tau_b_20_mid = integrals.tau_cluster_b(theta_s, phi_sub, dphi_sub,
                                               beta_s, beta_sub, phi_box, phi_b_20_mid,
                                               Du_dstream, Dv_stream,
                                               Du_ustream, Dv_ustream,
                                               lambda, box_height, gamma)
  if print_out_box then
    print('here_6_1')
  end
  local tau_b_03_mid = integrals.tau_cluster_b(theta_s, phi_sub, dphi_sub,
                                               beta_s, beta_sub, phi_box, phi_b_03_mid,
                                               Du_dstream, Dv_dstream, 
                                               Du_ustream, Dv_ustream,
                                               lambda, box_height, gamma)
  if print_out_box then
    print('here_6_1')
  end
  local tau_b_34_mid = integrals.tau_cluster_b(theta_s, phi_sub, dphi_sub,
                                               beta_s, beta_sub, phi_box, phi_b_34_mid,
                                               Du_dstream, Dv_stream,
                                               Du_ustream, Dv_ustream,
                                               lambda, box_height, gamma)
  if print_out_box then
    print('here_6_1')
  end

  local y_12 = integrals.geometry_b(tau_b_12_mid,phi_b_12)
  if print_out_box then
    print('here_6_1')
  end
  local x_20 = integrals.geometry_b(tau_b_20_mid,phi_b_20)
  if print_out_box then
    print('here_6_1')
  end
  local x_03 = integrals.geometry_b(tau_b_03_mid,phi_b_03)
  if print_out_box then
    print('here_6_1')
  end
  local y_34 = integrals.geometry_b(tau_b_34_mid,phi_b_34)
  if print_out_box then
    print('here_6_1')
  end

  local val = bernoulli_asym(tau_s_mid, eta_s_mid, froude) -- N-2 equations (solves at mid-points and alpha determined by fixing downstream location)

  -- Equation for linearised far field
  val[#val+1] = math.tan(lambda) - f_sq*lambda

  -- Ensure that each crest has zero velocity <- (#n_sub - 1) equations
  local k = 0
  for i = 1, (m-1) do
     k = k + n_sub[i]
     val[#val+1] = (theta_s_mid[k] + theta_s_mid[k+1])
  end

  if print_out_box then
    print('here_6_1')
  end
  -- Far field matching equations - 6 equations.
  local phi_lin_us = integrals.phi({beta_s[1], beta_s[2], (beta_s[1]+beta_s[2])/2.0}, phi_sub, dphi_sub)
  local phi_lin_ds = integrals.phi({beta_s[n], beta_s[n-1], (beta_s[n]+beta_s[n-1])/2.0}, phi_sub, dphi_sub)

  local u_lim_A = 1 - (Du_ustream*math.cos(lambda)*math.exp(-lambda*phi_lin_us[1]))
  local u_lim_B = 1 - (Du_ustream*math.cos(lambda)*math.exp(-lambda*phi_lin_us[2]))
  local v_lim_A = (Dv_ustream*math.sin(lambda)*math.exp(-lambda*phi_lin_us[1]))
  local v_lim_B = (Dv_ustream*math.sin(lambda)*math.exp(-lambda*phi_lin_us[2]))

  local theta_lim_A =  math.atan(v_lim_A / u_lim_A)
  local theta_lim_B =  math.atan(v_lim_B / u_lim_B)

  val[#val+1] = theta_s[1] - theta_lim_A
  val[#val+1] = theta_s[2] - theta_lim_B

  u_lim_A = 1 - (Du_dstream*math.cos(lambda)*math.exp(-lambda*phi_lin_ds[1]))
  u_lim_B = 1 - (Du_dstream*math.cos(lambda)*math.exp(-lambda*phi_lin_ds[2]))
  v_lim_A = -(Dv_dstream*math.sin(lambda)*math.exp(-lambda*phi_lin_ds[1]))
  v_lim_B = -(Dv_dstream*math.sin(lambda)*math.exp(-lambda*phi_lin_ds[2]))

  theta_lim_A =  math.atan(v_lim_A / u_lim_A)
  theta_lim_B =  math.atan(v_lim_B / u_lim_B)

  val[#val+1] = theta_s[n]   - theta_lim_A
  val[#val+1] = theta_s[n-1] - theta_lim_B

  local u_lim_C   = 1 - (Du_ustream*math.cos(lambda)*math.exp(-lambda*phi_lin_us[3]))
  local v_lim_C   = (Dv_ustream*math.sin(lambda)*math.exp(-lambda*phi_lin_us[3]))
  local tau_lim_C = math.log(u_lim_C * u_lim_C + v_lim_C * v_lim_C) / 2.0
  val[#val+1]     = tau_s_mid[1] - tau_lim_C

  u_lim_C     = 1 - (Du_dstream*math.cos(lambda)*math.exp(-lambda*phi_lin_ds[3]))
  v_lim_C     = -(Dv_dstream*math.sin(lambda)*math.exp(-lambda*phi_lin_ds[3]))
  tau_lim_C   = math.log(u_lim_C * u_lim_C + v_lim_C * v_lim_C) / 2.0
  val[#val+1] = tau_s_mid[n-1] - tau_lim_C
 
 
  -- Geometry equations - 3 equations.
  val[#val+1] = (x_20 + x_03) - box_width
  val[#val+1] = y_34 - math.abs(box_height)
  val[#val+1] = y_12 - math.abs(box_height)


  return val

end
