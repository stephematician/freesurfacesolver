--[[ file: pressure_residuals.lua

  Script containing the general functions required to do symmetric calculations
  for a pressure disturbance like e-x^2. based on box_residuals.lua as of
  5/2/2013

  Modified to conform to new module definitions (see rev03 of box
  problem). 21/5/2013

  Author        : Stephen Wade
  Date created  : 05/2/2013
  Last modified : 21/5/2013
]]--

require "surface"
require "integrals_pressure_homotopy"
require "grid_functions"

module("pressure_residuals", package.seeall)

bernoulli = function(u, v, eta, p, froude)
  local f_sq = froude * froude
  local val = {}
  local alpha =  (u[#u] * u[#u]) + (v[#v] * v[#v]) + (2.0 * p[#p] / f_sq)

  for i = 1, (#u)-1 do
    val[i] = (u[i] * u[i]) + (v[i] * v[i]) +
             (2.0 * (eta[i] + p[i]) / f_sq) -
             alpha
  end
  return val
end

compute_pressure = function(A, B, x)
  local sqrt_pi = math.sqrt(math.pi)
  local val = {}
  for i = 1, #x do
    val[i] = A * B * math.exp(-(B * B) * x[i] * x[i]) / sqrt_pi
  end

  return val
end

compute_buvxe_hs = function(unpacked)

  -- v on free surface
  local v_s    = unpacked.V_S

  -- Parameters of farstream
  local lambda     = unpacked.LAMBDA
  local D          = unpacked.D

  -- Parameters of grid
  local n_sub      = unpacked.N_SUB
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local homotopy_s = unpacked.HOMOTOPY_S

  local n = #v_s

  local beta_s =  grid_functions.calculate_beta(n_sub)

  local u_s_mid = integrals_pressure.u_symmetric_s(v_s, phi_sub, dphi_sub,
                                                   beta_s, beta_sub, D,
                                                   lambda, homotopy_s)
  
  local v_s_mid = {}
  for i = 1, n-1 do
     v_s_mid[i] = (v_s[i] + v_s[i+1]) / 2.0
  end

  local x_s_mid, 
          x_lhs = integrals_pressure.x_symmetric_s(u_s_mid, v_s_mid,
                                                   beta_s, beta_sub, phi_sub, dphi_sub, 
                                                   homotopy_s)

  local eta_s_mid,
          eta_lhs = integrals_pressure.y_symmetric_s(u_s_mid, v_s_mid,
                                                     beta_s, beta_sub, phi_sub, dphi_sub,
                                                     homotopy_s)
 
  for j = 1, (n-1) do
    eta_s_mid[j] = (eta_s_mid[j] - eta_s_mid[#eta_s_mid])
  end

  return beta_s, u_s_mid, v_s_mid, x_s_mid, x_lhs, eta_s_mid, eta_lhs

end

compute_output_hs = function(unpacked)
  local t = {}

  local t = {}
    
  local beta_s, u_s_mid, v_s_mid,
                x_s_mid,   x_lhs,
              eta_s_mid, eta_lhs = compute_buvxe_hs(unpacked)

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
  local phi_s_mid  = integrals_pressure.phi(beta_s_mid, beta_sub, phi_sub, dphi_sub, homotopy_s)

  -- Free surface variables
  t.u_s_mid = u_s_mid
  t.v_s_mid = v_s_mid
  t.x_s_mid = x_s_mid
  t.eta_s_mid = eta_s_mid
  t.phi_s_mid = phi_s_mid

  -- Parameters of flow
  t.A = unpacked.A
  t.B = unpacked.B
  t.F = unpacked.FROUDE

  -- Parameters of far stream
  t.D     = unpacked.D
  t.lambda = unpacked.LAMBDA

  -- Parameters of grid
  t.n_sub       = unpacked.N_SUB
  t.phi_sub     = unpacked.PHI_SUB
  t.dphi_sub    = unpacked.DPHI_SUB
  t.homotopy_s  = unpacked.HOMOTOPY_S

  -- Now required 23/06/2014 for sensitivity analysis
  t.cluster_val = unpacked.CLUSTER_VAL
  
  return t
  
end

unit_clamp = function(g)
  return math.max(math.min(g,1), 0)
end

residual_hs = function(unknowns, extra)
  local fs = extra[1]
  
  local unpacked = surface.unpack_vec(unknowns,
                                    fs.initial,
                                    fs.is_free,
                               fs.is_continued,
                                  fs.ind_table)

  -- v on the free surface
  local v_s    = unpacked.V_S

  -- Flow parameters
  local froude     = unpacked.FROUDE
  local eta0       = unpacked.ETA0

  -- Pressure disturbance parameters
  local A        = unpacked.A
  local B        = unpacked.B

  -- Far-stream parameters
  local lambda     = unpacked.LAMBDA
  local D          = unpacked.D

  -- Grid parameters
  local beta_sub   = unpacked.BETA_SUB
  local phi_sub    = unpacked.PHI_SUB
  local dphi_sub   = unpacked.DPHI_SUB
  local n_sub      = unpacked.N_SUB
  local homotopy_s = unpacked.HOMOTOPY_S
  
  -- Optional flow parameters
  local u_min       = unpacked.U_MIN
  local cluster_val = unpacked.CLUSTER_VAL

--  local F_SQ = FROUDE * FROUDE

  if print_out then
    print('here_6_2')
  end

  local beta_s, u_s_mid, v_s_mid,
                x_s_mid,   x_lhs,
              eta_s_mid, eta_lhs = compute_buvxe_hs(unpacked)

  local beta_s_mid = {}
  for i = 1, (#beta_s)-1 do
    beta_s_mid[i] = (beta_s[i] + beta_s[i+1]) / 2.0
  end

  if print_out then
    print('here_6_3')
  end

  local phi_s_mid  = integrals_pressure.phi(beta_s_mid, beta_sub, phi_sub, dphi_sub, homotopy_s)

  if print_out then
    for kk,vv in pairs(phi_sub) do
      print(kk,vv)
    end
    for kk,vv in pairs(dphi_sub) do
      print(kk,vv)
    end
  end

  local n = (#v_s)

  local p_s = compute_pressure(A, B, x_s_mid)
  
  --local p_0 = p[n-1]

  --local u_0 = u_mid[n-1]
  --local v_0 = v_mid[n-1]

  local val = bernoulli(u_s_mid, v_s_mid, eta_s_mid, p_s, froude) --, alpha) -- N-2 equations

  -- Equations for crest(s)
  local k = 0
  for i = 1, (#n_sub-1) do
     k = k + n_sub[i]
     val[#val+1] = (v_s_mid[k] + v_s_mid[k+1])
  end

  -- Far field matching equations
  -- Same scheme as box integral method
  -- (2 based on the 'unknown' - tau/v, and 1 based on the 'known', theta/u)
  phi_s_ff_a  = phi_s_mid[n-1]
  -- u_s_ff_a    = 1 - (D*math.cos(lambda)*math.exp(-lambda*phi_s_ff_a))
  v_s_ff_a    =    -(D*math.sin(lambda)*math.exp(-lambda*phi_s_ff_a))
  -- val[#val+1] = u_s_mid[n-1] - u_s_ff_a
  val[#val+1] = v_s_mid[n-1] - v_s_ff_a

  phi_s_ff_b  = phi_s_mid[n-2]
  v_s_ff_b    = -(D*math.sin(lambda)*math.exp(-lambda*phi_s_ff_b))
  val[#val+1]  = v_s_mid[n-2] - v_s_ff_b

  -- Equation for linearised far field
  val[#val+1] = math.tan(lambda) - froude*froude*lambda

  -- Geometry equations
  val[#val+1] = eta_s_mid[1] - eta0

  -- Optional parameters
  if u_min then
    local integral_u_s_min = 1/0
    for kk, vv in pairs(u_s_mid) do
      if vv < integral_u_s_min then
        integral_u_s_min = vv
      end
    end
    val[#val+1] = u_min - integral_u_s_min
  end

  if cluster_val then
    local k = 1
    for i = 1, (#n_sub) do
      local u_k = unit_clamp(u_s_mid[k])
      if print_out then print(u_k) end
      val[#val+1] = dphi_sub[i] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n) -- 129
      k = k + n_sub[i]
    end
    assert(k == n, 'Failed assertion k == n for some thing')
    local u_k = unit_clamp(u_s_mid[#u_s_mid])
    if print_out then print(u_k) end
    val[#val+1] = dphi_sub[#dphi_sub] - (((cluster_val * math.pow(1-u_k,2)) + (phi_sub[#phi_sub] * math.pow(u_k,2)))/n) -- 130
  end

  if print_out then
    for kk, vv in pairs(val) do
      print(kk, vv)
    end
    io.flush()
  end


  return val
end