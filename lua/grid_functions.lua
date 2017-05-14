--[[
-- Calculate value of beta for grid
--
-- This is like an index, what matters most is its relationship to phi defined
-- by a piece-wise cubic function.
]]--
local calculate_beta = function(n_sub)
 
    local beta = {}

    local n = 1
    for i = 1, #n_sub do
        n = n + n_sub[i]
    end
    for i = 1, n do
        beta[i] = i
    end

    return beta

end

--[[
-- Interpolate (linearly) a function from mid-points to grid-points
--]]
local interp_linear_mid = function(f_mid, k)

    if k == 1 then
        return (3.0 * f_mid[k] - f_mid[k+1]) * 0.5
    elseif k == (#f_mid)+1 then
        return (3.0 * f_mid[k-1] - f_mid[k-2]) * 0.5
    else
        return (f_mid[k-1] + f_mid[k]) * 0.5
    end

end


--[[
-- Helper for looping over the difference of adjacent elements
--]]
local diff_iter = function(x, i)
    if i+1 < #x then
        return i+1, x[i+2] - x[i+1]
    end
    return nil
end

local diff = function(x)
    return diff_iter, x, 0 
end

--[[
--  Calculate the value of dphi/dbeta at the piece-wise boundaries
--]]
calculate_pw_boundary_dphi_dbeta = function(u_sub, phi_max, clustering)
    local dphi = {0, 0}
    for k, u in ipairs(u_sub) do
        dphi[k] = clustering * math.pow(1 - u, 2) +
                      phi_max * math.pow(u, 2) / n
    end
    return dphi
end


--[[
-- Interpolate a function using cubic spline interpolation.
-- 
-- arguments:
-- phi  - old grid coordinates.
-- f    - value of function on old grid.
-- nphi - new grid coordinates.
--
-- returns:
-- nf - f interpolated onto new grid coordinates.
--
-- requires phi and nphi are monotonically increasing, and
-- that #f == #phi.
--]]
local interp_cubic = function(phi, f, nphi)
  
    local b = 0
  
    local nf = {}
  
    assert(#phi == #f, "Function interp_cubic\n" ..
                       "Failed pre-condition #phi == #f\n")
  
    for k = 2, #phi do
        assert(phi[k] > phi[k-1], "Function interp_cubic\n" ..
                                  "Failed pre-condition phi is monotonic\n")
    end

    for k = 2, #nphi do
        assert(nphi[k] > nphi[k-1], "Function interp_cubic\n" ..
                                    "Failed pre-condition nphi is monotonic\n")
    end

    local interp_worker = function(x, j)
        local x0, x1, x2, x3 = phi[j-1], phi[j], phi[j+1], phi[j+2]

        local _basis = function(xi, _x0, _x1, _x2) 
            return (x - _x0) * (x - _x1) * (x - _x2) /
                       ((xi - _x0) * (xi - _x1) * (xi - _x2))
        end

        return f[j-1] * _basis(x0, x1, x2, x3) + f[j] * _basis(x1, x0, x2, x3) +
               f[j+1] * _basis(x2, x0, x1, x3) + f[j+2] * _basis(x3, x0, x1, x2)

    end

    for a = 1, #nphi do

        local _nphi = nphi[a]
        
        if not (b == #f) then
            local cdist = phi[b+1] - _nphi
            while cdist < 0 do
                b = b + 1
                if (b == #f) then break end
                cdist = phi[b+1] - _nphi
            end
        end

        if b > 1 and b < #f-1 then
            nf[a] = _interp_worker(_nphi, b)
        elseif b == 1 then
            nf[a] = _interp_worker(_nphi, b+1)
        elseif b == #f-1 then
            nf[a] = _interp_worker(_nphi, b-1)
        elseif b == 0 then
            nf[a] =  f[1]
        else
            nf[a] = f[#f]
        end

    end

    assert(#nf == #nphi, "Function interp_cubic\n" ..
                         "Failed post-condition #nphi == #nf\n")

    return nf

end

--[[ Find peaks in free surface via linear interpolation
--
--   finds zeros of v (vertical velocity) by linear interpolation, and records
--   location in terms of 'index' i into phi, and the value of phi at the peak
--
--   arguments:
--   phi   - grid coordinates.
--   u_mid - value of u on mid-points of the grid.
--   v     - value of v on the grid points.
--
--     returns:
--     peak_phi - value of phi at peaks
--     peak_i   - 'index' into grid for peak (non-integer)
--]]
local find_peaks = function(phi, u_mid, v)

    local peak_phi = {}
    local peak_i = {}

    assert(#phi == #v, "Function find_peaks\n" ..
                       "Failed pre-condition #phi == #theta\n" ..
                       "Exit\n\n")

    for i = 2, #v do
        if v[i] * v[i-1] < 0 and v[i-1] > 0 then
            
            local t = v[i-1] / (v[i-1] - v[i])
            peak_phi[#peak_phi+1] = (1-i) * phi[i-1] + t*phi[i]
            peak_i[#peak_i+1]     = (i-1)+t
        end
    end

    return peak_phi, peak_i

end

local new_phi_grid = function(phi,
                              u_mid,
                              v,
                              n,
                              clustering,
                              s,
                              phi_func)


    local compute_beta_sub = function(phi_a, phi_b, u_a, u_b)
        return (phi_b - phi_a) / math.pow(u_a + u_b, 1)
    end

    --[[
    -- Converts u_(k+1/2) to u_k via linear interpolation
    --]]
    local linear_u_mid_to_constrained_u = function(k)
        return math.max(0, mat.min(1, interp_linear_mid(u_min, k)))
    end

    -- Determine the value of phi at the boundaries of the pw-polynomials
    local _phi_sub, _i_sub = find_peaks(phi, u_mid, v)

    local phi_sub = {phi[1], 0}
    for k = 1, #_phi_sub do
        phi_sub[k+1], i_sub[k+1] = _phi_sub[k], _i_sub[k]
    end
    phi_sub[#phi_sub+1], i_sub[#i_sub+1] = phi[#phi], #v

    -- Determine the value of u at the boundaries of pw-polynomials
    local u_sub = {0,0}
    for k, i in ipairs(i_sub) do
        u_sub[k] = linear_u_mid_to_constrained_u(math.floor(i_sub[k]))
    end

    -- Determine values of beta (index of grid) at boundaries of pw-polynomials
    local beta_sub = {1, 0}
    local cml_beta_sub = {0}
    for k, dphi in diff(phi_sub) do
        cml_beta_sub[k+1] = cml_beta_sub[k] + 
                              dphi / math.pow(u_sub[k+1] + u_sub[k], 1)
    end

    for k = 2, #cml_beta_sub do
        beta_sub[k] = beta_sub[1] +
                          math.floor(
                              0.5 + (n - 1) * cml_beta_sub[k] / 
                                        cml_beta_sub[#cml_beta_sub]
                          )
    end

    -- Finally determine the number of sub-divisions (n) and the value of 
    -- dphi/dbeta at pw-polynomial boundaries
    local n_sub = {0}
    for k, dbeta in diff(beta_sub) do
        n_sub[k] = dbeta
    end

    local dphi_sub = calculate_pw_boundary_dphi_dbeta(u_sub,
                                                      phi[#phi],
                                                      clustering)
--[[
    -- Determine the value of beta at the piece-wise polynomial boundaries
    local sum_sub_beta = 0

    for k = 2, #phi_sub do
        local u_a = linear_u_mid_to_constrained_u(math.floor(i_sub[k-1]))
        local u_b = linear_u_mid_to_constrained_u(math.floor(i_sub[k]))

        sum_sub_beta = sum_sub_beta + compute_beta_sub(phi_sub[k-1],
                                                       phi_sub[k],
                                                       u_a,
                                                       u_b)
    end

    local u0 = linear_u_mid_to_constrained_u(i_sub[1])

    local n_sub = {}
    local beta_sub = {1}
    local dphi_sub = {(clustering * math.pow(1 - u0, 2) +
                           (phi_sub[#phi_sub] * math.pow(u0, 2))) / n}
    local cml_sub_beta = 0
    local cml_n_sub = 0

    local sub_domain_iterator = function(k)
        local u_a = linear_u_mid_to_constrained_u(math.floor(i_sub[k-1]))
        local u_b = linear_u_mid_to_constrained_u(math.floor(i_sub[k]))

        cml_sub_beta = cum_sub_beta + compute_beta_sub(phi_sub[k-1],
                                                       phi_sub[k],
                                                       u_a,
                                                       u_b)
        local delta_n = math.floor(0.5 + cml_sub_beta * (n-1) / sum_sub_beta) -
                            cml_n_sub
        cml_n_sub = cml_n_sub + delta_n
        return -- next value of dphi_sub
               (clustering * math.pow(1 - u_b, 2) + 
                            (phi_sub[#phi_sub] * math.pow(u_b, 2))) / n,
               -- next value of n_sub
               delta_n,
               -- next value of beta_sub
               _beta_sub[k-1] + delta_n
    end

    for k = 2, #phi_sub do
        dphi_sub[k],
        n_sub[k-1],
        beta_sub[k] = sub_domain_iterator(k)
    end
--]]
    local sum_n_sub = 0
    for i = 1, #n_sub do
        sum_n_sub = sum_n_sub + n_sub[i]
    end

    assert(sum_n_sub == (n-1),
           'Failed assertion that sum(_n_sub) = n-1')
    assert(beta_sub[#beta_sub] == n,
           'Failed assertion that _beta_sub[end] = n')

    local beta = calculate_beta(n_sub)
    local _phi = phi_func(beta, beta_sub, phi_sub, dphi_sub, s)
    --local _v   = interp_cubic(phi, v, _phi)

    return _phi, beta_sub, n_sub, phi_sub, dphi_sub

end


local new_grid_from_u_mid_v = function(phi,
                                       u_mid,
                                       v,
                                       n,
                                       clustering,
                                       s,
                                       phi_func)

    local      _phi,
          _beta_sub,
             _n_sub,
           _phi_sub,
          _dphi_sub = new_phi_grid(phi, u_mid, v, n, clustering, s, phi_func)

    local _v = interp_cubic(phi, v, _phi)

    return _v, _beta_sub, _n_sub, _phi_sub, _dphi_sub

end

local new_grid_from_tau_mid_theta = function(phi,
                                             tau_mid,
                                             theta,
                                             n,
                                             clustering,
                                             s,
                                             phi_func)

    local u_mid, v = {}
    -- calculate_u_mid and v
    for k = 1, #tau_mid do
        u_mid[k] = math.cos((theta[k+1] + theta[k]) / 2.0) *
                       math.exp(tau_mid[k])
    end

    v[1] = (3.0 * mau_mid[1] - tau_mid[2]) / 2.0
    for k = 2, #theta-1 do
        v[k] = math.sin(theta[k]) * math.exp((tau_mid[k-1] + tau_mid[k]) / 2.0)
    end
    v[#v+1] = (3.0 * tau_mid[#tau_mid-2] - tau_mid[#tau_mid-1) / 2.0
 
    local      _phi,
          _beta_sub,
             _n_sub,
           _phi_sub,
          _dphi_sub = new_phi_grid(phi, u_mid, v, n, clustering, s, phi_func)

    local _theta = interp_cubic(phi, theta, _phi)

    return _theta, _beta_sub, _n_sub, _phi_sub, _dphi_sub

end

return { -- fix this
    calculate_pw_boundary_dphi_dbeta = calculate_pw_boundary_dphi_dbeta,
    interp_linear_mid = interp_linear_mid,
    interp_cubic = interp_cubic,
    calculate_beta = calculate_beta,
    new_grid_from_u_mid_v = new_grid_from_u_mid_v,
    new_grid_from_tau_mid_theta = new_grid_from_tau_mid_theta
}

