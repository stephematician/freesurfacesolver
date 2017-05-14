local grid = require "grid_functions"
local c_integrals = require "integrals_sbox_homotopy"

local tau_s_mid_symmetric = function(unpacked, beta_s)

    local _beta_s = beta_s 
    if not _beta_s then
        _beta_s = grid.calculate_beta(unpacked.N_SUB)
    end


    return c_integrals.tau_smoothbox_symmetric_s(unpacked.THETA_S,
                                                 unpacked.PHI_SUB,
                                                 unpacked.DPHI_SUB,
                                                 _beta_s
                                                 unpacked.BETA_SUB,
                                                 unpacked.A,
                                                 unpacked.B,
                                                 unpacked.L,
                                                 unpacked.PHI_C,
                                                 unpacked.D,
                                                 unpacked.LAMBDA,
                                                 unpacked.GAMMA,
                                                 unpacked.HOMOTOPY_S)
end

local theta_s_mid_symmetric = function(unpacked)

    local theta_s_mid = {0}
    for i = 1, n-1 do
         theta_s_mid[i] = (unpacked.THETA_S[i] + unpacked.THETA_S[i+1]) / 2.0
    end

    return theta_s_mid

end

local y_s_mid = function(unpacked, beta_s, theta_s_mid, tau_s_mid)

    local _beta_s = beta_s
    local _theta_s_mid = theta_s_mid
    local _tau_s_mid = tau_s_mid

    if not _beta_s then
        _beta_s = beta_s = grid.calculate_beta(unpacked.N_SUB)
        _theta_s_mid = theta_s_mid_symmetric(unpacked)
        _tau_s_mid = tau_s_mid_symmetric(unpacked, _beta_s)
    end

    return c_integrals.y_s(       _theta_s_mid,        _tau_s_mid,
                                       _beta_s, unpacked.BETA_SUB,
                              unpacked.PHI_SUB, unpacked.DPHI_SUB,
                           unpacked.HOMOTOPY_S)

end

local x_s_mid = function(unpacked, beta_s, theta_s_mid, tau_s_mid)

    local _beta_s = beta_s
    local _theta_s_mid = theta_s_mid
    local _tau_s_mid = tau_s_mid

    if not _beta_s then
        _beta_s = beta_s = grid.calculate_beta(unpacked.N_SUB)
        _theta_s_mid = theta_s_mid_symmetric(unpacked)
        _tau_s_mid = tau_s_mid_symmetric(unpacked, _beta_s)
    end

    return c_integrals.x_s(       _theta_s_mid,        _tau_s_mid,
                                       _beta_s, unpacked.BETA_SUB,
                              unpacked.PHI_SUB, unpacked.DPHI_SUB,
                           unpacked.HOMOTOPY_S)

end

local phi_s_mid = function(unpacked, beta_s_mid)

    local _beta_s_mid = beta_s_mid
    if not _beta_s_mid then
        local beta_s = grid.calculate_beta(unpacked.N_SUB)
        for k = 1, #beta_s-1 do
            _beta_s_mid[k] = 0.5 * (beta_s[k] + beta_s[k+1])
        end
    end

    return c_integrals.phi(_beta_s_mid,
                           unpacked.BETA_SUB,
                           unpacked.PHI_SUB,
                           unpacked.DPHI_SUB,
                           unpacked.HOMOTOPY_S)

end

return {
    tau_s_mid_symmetric = tau_s_mid_symmetric,
    theta_s_mid = theta_s_mid_symmetric,
    x_s_mid = x_s_mid,
    y_s_mid = y_s_mid,
    phi_s_mid = phi_s_mid
}

