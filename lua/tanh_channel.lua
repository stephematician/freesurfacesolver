local generic_regrid = require "topography.regrid"
local generic_load   = require "topography.load"
local integrals = require "integrals_sbox_homotopy"
local interface = require "tanh_channel.integral_interface"
local grid = require "grid_functions"

local residual_symmetric
local output_symmetric


local regrid_symmetric = function(_surf, _prefix)

    return generic_regrid.regrid(_surf,
                                 integrals.phi_s,
                                 interface.tau_s_mid_symmetric,
                                 _prefix)

end


local load_last_symmetric = function(file,
                                     init_tables,
                                     n,
                                     output_data,
                                     _prefix)

    return generic_load.load_last(file,
                                  integrals.phi_s,
                                  init_tables,
                                  n,
                                  residual_symmetric,
                                  output_symmetric,
                                  output_data,
                                  _prefix)

end


local unit_clamp = function(g)
    return math.max(0,math.min(1,g))
end


residuals_symmetric = function(unknowns, expand)

    local unpacked = expand(unknowns)

    local beta_s = grid.calculate_beta(unpacked.N_SUB)
    local tau_s_mid = interface.tau_s_mid_symmetric(unpacked, beta_s)
    local theta_s_mid = interface.theta_s_mid_symmetric(unpacked)

    local eta_s_mid, eta_lhs = interface.y_s_mid(unpacked,
                                                 beta_s,
                                                 theta_s_mid,
                                                 tau_s_mid)

    for k = 1, #eta_s_mid do
        eta_s_mid[k] = eta_s_mid[k] - eta_s_mid[#eta_s_mid]
    end

    --[[
    --  N-2 equations from bernoulli's condition on the free-surface
    --]]
    local val = bernoulli(tau_s_mid, eta_s_mid, unpacked.FROUDE)

    --[[
    --  #crests equations - for aligning piece-wise polynomial sub-domains with 
    --                      crests
    --]] 
    local k_crest = 0
    for k = 1, (#n_sub-1) do
        k_crest = k_crest + n_sub[k]
        val[#val+1] = (theta_s_mid[k_crest] + theta_s_mid[k_crest+1])
    end
    
    
    local beta_s_mid = {0, 0}
    for k = 1, #beta_s-1 do
        beta_s_mid[k] = 0.5 * (beta_s[k] + beta_s[k+1])
    end

    local n = #theta_s_mid+1

    local phi_s_ff  = interface.phi_s_mid_symmetric(unpacked,
                                                    {beta_s_mid[n-1],
                                                     beta_s_mid[n-2]})

    --[[
    --  2 far-field matching equations - to ensure linear system matches at final
    --                                   two grid points
    --]]
    local decay_ff_1 = math.exp(-unpacked.LAMBDA*phi_s_ff[1])
    val[#val+1] = math.atan(unpacked.D * math.sin(unpacked.LAMBDA) * 
                                decay_ff_1 / 
                                (1 + unpacked.D * math.cos(unpacked.LAMBDA) *
                                         decay_ff_1)) - theta_s_mid[n-1]
    local decay_ff_2 = math.exp(-unpacked.LAMBDA*phi_s_ff[2])
    val[#val+1] = math.atan(unpacked.D * math.sin(unpacked.LAMBDA) * 
                                decay_ff_2 / 
                                (1 + unpacked.D * math.cos(unpacked.LAMBDA) *
                                         decay_ff_2)) - theta_s_mid[n-2]
    --[[
    --  Linearised far-field condition
    --]]
    val[#val+1] = math.tan(unpacked.LAMBDA) -
                      math.pow(unpacked.FROUDE, 2.0) * unpacked.LAMBDA

    --[[
    --  Geometry equations
    --]]
    local x_b = interface.x_b_symmetric(unpacked,
                                        {unpacked.PHI_C},
                                        beta_s)

    --[[
    -- 1 equation to ensure topography width is correct
    --]]
    val[#val+1] = x_b[1] - (unpacked.L * 0.5)

    --[[
    -- 1 equation to define eta0 parameter
    --]]
    val[#val+1] = eta_s_mid[1] - unpacked.ETA0

    --[[
    -- #n_sub + 1 equations to determine dphi at piecewise polynomial end-points
    --]]
    if unpacked.CLUSTER_VAL then
        local u_mid = {0, 0}
        for k, j in ipairs(unpacked.BETA_SUB) do
            if j < #theta_s_mid then
                u_mid[j], u_mid[j+1] = unit_clamp(math.cos(theta_s_mid[j]) * 
                                                      math.exp(tau_s_mid[j])),
                                       unit_clamp(math.cos(theta_s_mid[j+1]) *
                                                      math.exp(tau_s_mid[j+1]))
            elseif j == #theta_s_mid then
                u_mid[j-1], u_mid[j] = unit_clamp(math.cos(theta_s_mid[j-1]) * 
                                                      math.exp(tau_s_mid[j-1])),
                                       unit_clamp(math.cos(theta_s_mid[j]) *
                                                      math.exp(tau_s_mid[j]))
            else
                error([=[Unexpected index found in grid specification]=] ..
                      [=[beta_sub]=])
            end
            u_sub[k] = grid.interp_linear_mid(u_mid, j)
        end
        local dphi_sub = calculate_pw_boundary_dphi_dbeta(
                             u_sub,
                             unpacked.PHI_SUB[#unpacked.PHI_SUB],
                             unpacked.CLUSTER_VAL
                         )
        for k, dphi in ipairs(unpacked.DPHI_SUB) do
            val[#val+1] = dphi - dphi_sub[k]
        end
    end

    --[[
    -- 1 equation (optional) if u_min parameter exists
    --]]
    if unpacked.U_MIN then
        local u_min = 1/0
        for k, tau in ipairs(tau_s_mid) do
            u_min = math.min(u_min, math.cos(theta_s_mid[k]) * math.exp(tau))
        end
        val[#val+1] = unpacked.U_MIN - u_min
    end

end


output_symmetric = function()
end

return {
    regrid_symmetric = regrid_symmetric,
    load_last_symmetric = load_last_symmetric,
    residual_symmetric = residual_symmetric,
    output_symmetric = output_symmetric
}
