local grid = require "grid_functions"

local DEBUG_OUTPUT = function(_prefix,
                              ...)
    
    local prefix = _prefix
    if VERBOSE then
        if not prefix then
            prefix = ''
        end
        io.write(prefix, ...)
    end
    io.flush()

end

--[[
-- Re-grid worker for topography problems
--
-- Calculates the form of the grid for a surface, and if necessary, constructs
-- a new surface on the new grid when required. This new surface is then
-- checked for the direction that the continuation algorithm proceeds in.
-- 
-- params:
--  - _surf: the original surface
--  - phi_func: a function which calculates phi from beta (and the grid
--              beta_sub, phi_sub, dphi_sub and clustering parameters from the
--              original surface)
--  - tau_mid_func: a function which calculates tau_mid on the free surface 
--                  from the surface specification data. 
--  - prefix: a prefix to attach to each line out verbose output
--
-- returns:
--  - the original or new surface
--  - the direction of continuation
--]]
local regrid = function(_surf,
                        phi_func,
                        tau_mid_func,
                        _prefix)

    local prefix = _prefix
    if not prefix then
        prefix = ''
    end

    local _initial      = _surf.expand()
    local _is_free      = _surf.get_is_free()
    local _is_continued = _surf.get_is_continued()

    local _n_sub = _initial.N_SUB

    local _beta_s    = grid.calculate_beta(_initial.N_SUB)
    local _tau_s_mid = tau_mid_func(_initial, _beta_s)
    local _phi_s     = phi_func(_beta_s,
                                _initial.BETA_SUB,
                                _initial.PHI_SUB,
                                _initial.DPHI_SUB,
                                _initial.CLUSTER_VAL)
    DEBUG_OUTPUT(prefix,
                 [=[    regrid()]=], '\n')
    DEBUG_OUTPUT(prefix,
                 [=[        Calculating grid for current surface.]=],
                 '\n')
 
    _initial.THETA_S,
    _initial.BETA_SUB,
    _initial.N_SUB,
    _initial.PHI_SUB,
    _initial.DPHI_SUB = grid.new_grid_from_tau_mid_theta(
                            _phi_s,
                            _tau_s_mid
                            _initial.THETA_S,
                            #(_initial.THETA_S),
                            _initial.CLUSTER_VAL,
                            _initial.HOMOTOPY_S,
                            phi_func
                        )

    -- Determine if we need to re-grid
    local need_regrid = false

    if #(_initial.N_SUB) == _n_sub then
        for k = 1, #(_initial.N_SUB) do
            if not (_initial.N_SUB[k] == _n_sub[k]) then
                need_regrid = true
                break
            end
        end
    else
        need_regrid = true
    end

    if need_regrid then
        if VERBOSE then
            io.write(prefix,
                     [=[    regrid()]=], '\n')
            io.write(prefix,
                     [=[        Now require new surface continuation object.]=],
                     '\n')
            io.flush()
        end
    
        _is_free.BETA_SUB, _is_continued.BETA_SUB = {}, {}
        for k = 1, #_initial.BETA_SUB do
            _is_free.BETA_SUB, _is_continued.BETA_SUB = false, false
        end

        _is_free.N_SUB, _is_continued.N_SUB = {}, {}
        for k = 1, #_initial.N_SUB do
            _is_free.N_SUB[k], _is_continued.N_SUB[k] = false, false
        end

         _is_free.PHI_SUB,  _is_continued.PHI_SUB = {false}, {false}
        _is_free.DPHI_SUB, _is_continued.DPHI_SUB =  {true}, {false}
        for k = 2, #_initial.PHI_SUB do
             _is_free.PHI_SUB[k],  _is_continued.PHI_SUB[k] = true, false
            _is_free.DPHI_SUB[k], _is_continued.DPHI_SUB[k] = true, false
        end
                    _is_free.PHI_SUB[#(_is_free.PHI_SUB)+1],
          _is_continued.PHI_SUB[#(_is_continued.PHI_SUB)+1] = false, false
                  _is_free.DPHI_SUB[#(_is_free.DPHI_SUB)+1],
        _is_continued.DPHI_SUB[#(_is_continued.DPHI_SUB)+1] =  true, false

        if VERBOSE then
            io.write(prefix,
                     [=[    regrid()]=], '\n')
            io.write(prefix,
                     [=[        Calling surface copy constructor.]=],
                     '\n')
            io.flush()
        end

        local regrid_surf = _surf.copy(_initial,
                                       _is_free,
                                       _is_continued)

        if VERBOSE then
            io.write(prefix,
                     [=[    regrid()]=], '\n')
            io.write(prefix,
                     [=[        Solving surface on new grid with Newtons.]=],
                     '\n')
            io.flush()
        end
 
        local result = regrid_surf.solve_and_reset()

        assert(result == continuator.CONVERGED,
               [=[topography_grid() could not solve for new grid]=])

               if VERBOSE then
            io.write(prefix,
                     [=[    regrid()]=], '\n')
            io.write(prefix,
                     [=[        Performing continuation step to establish ]=],
                     [=[direction]=], '\n')
            io.flush()
        end
 
        local surface_data_a = regrid_surf.expand()
        local _direction = _surf.get_direction()
        regrid_surf.set_step(regrid_surf.get_step() * 1E-2)
        
        result = regrid_surf.progress_interp(0, -- dummy value
                                             _direction)
        assert((result == continuator.CONVERGED) or 
                    (result == continuator.BIFURCATION),
                [=[topography_regrid() could not perform continuation step.]=])
 
        local surface_data_b = regrid_surf.expand()

        local switched = _surf.test_switch(surface_data_a, surface_data_b)

        if switched then

            if VERBOSE then
                io.write(prefix,
                         [=[    regrid()]=], '\n')
                io.write(prefix,
                         [=[        Re-grid has caused direction change, ]=],
                         [=[will try to perform backwards step.]=]'\n')
                io.flush()
            end
            if _direction == continuator.FORWARD then
                _direction = continuator.BACKWARD
            else
                _direction = continuator.FORWARD
            end

            regrid_surf.flush_history()
            regrid_surf.set_step(regrid_surf.get_step() * 1E-2)
            
            result = regrid_surf.progress_interp(0, -- dummy value
                                                 _direction)
            assert((result == continuator.CONVERGED) or 
                        (result == continuator.BIFURCATION),
                    [=[topography_regrid() could not perform continuation ]=],
                    [=[step.]=])
        end

        return regrid_surf, regrid_surf.get_direction()

    else

        return _surf, _surf.get_direction()

    end

end

return {
    regrid = regrid   
}

