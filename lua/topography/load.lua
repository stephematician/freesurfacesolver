local surface = require "surface"
local grid = require "grid_functions"

local DEBUG_OUTPUT = function(_prefix, ...)

    local prefix = _prefix
    
    if VERBOSE then
        if not prefix then
            prefix = ''
        end
        io.write(prefix, ...)
    end
    io.flush()

end

local last_surface_from_matlab = function(file)

    local of = assert(io.open(file, "r"))

    -- First step, find %% 000... 
    of:seek("set", string.len("%% "))
    -- The last surface in the file
    local n = of:read("*n")

    -- Parse every line that is of the form
    -- o.*{n} = [ .... ];
    -- o.*{n} = yyyy;

    local surf = {}

    local valid_line_pattern = 'o%.(.+){' .. n .. '}'
    local value_pattern = '%-?%d+%.?%d*e?-?%d*'

    local line = of:read("*line")

    while (not (line == nil)) do

        local key_start, key_end = line:find(valid_line_pattern)
        
        if key_start then
            local key = line:sub(key_start, key_end)
            
            local value_start, value_end = line:find('%[', key_end + 1)

            local i = 1
            local value = nil

            if not (value_start == nil) then
                -- we have a table
                surf[key] = {}
                repeat
                    value_start, value_end = line:find(value_pattern,
                                                       value_end + 1)
                    surf[key][count] = tonumber(line:sub(value_start,
                                                         value_end))
                    i = i + 1

                until value_end == nil
            else
                -- we should have a single number
                value_start, value_end = line:find(value_pattern, key_end + 1)
                surf[key] = tonumber(line:sub(value_start, value_end))
            end
        end

        io.flush()
        line = of:read("*line")

    end

    of:close()

    return surf

end

--[[
-- Loads the last surface from a specified file, interpolates it onto a
-- desired grid size, then solves on the new grid and returns the surface
--
-- arguments
--
--    - file; name of the file to load the last surface from
--    - init_tables; function returning tables for initial, free variable and 
--        continuation variable specification for a surface with two parameters,
--        the grid size and number of crests (?)
--    - n; the desired grid size for the return surface
--
-- returns
--
--    -
--
--]]
local load_last = function(file,
                           phi_func,
                           init_tables,
                           n,
                           surface_residual,
                           surface_output,
                           output_data,
                           _prefix)

    local prefix = _prefix
    if not prefix then
        prefix = ''
    end

    local last_surface = last_surface_from_matlab(file)

    -- Use the init_tables function to create a new dummy surface description
    -- for the desired grid size
    --

    assert(not(last_surface == nil),
           [=[load_last() could not load final surface from specified file]=])

    local _n = #(last_surface.THETA_S_MID) + 1
    local _beta_s = grid.calculate_beta(last_surface.N_SUB)
    local _phi_s  = phi_func(_beta_s,
                             last_surface.BETA_SUB,
                             last_surface.PHI_SUB,
                             last_surface.DPHI_SUB,
                             last_surface.HOMOTOPY_S)

    local new_initital, _, _ = init_tables_fn(_n, #(last_surface.N_SUB))

    for k, v in pairs(new_initial) do
        if type(v) == "number" then
            assert(
                not (last_surface[k] == nil),
                [=[load_last() could not find data in specified file for ]=] .. k
            )
            new_initial[k] = last_surface[k]
        else
            assert(type(v) == "table",
                   [=[load_last() initial table must not contain ]=] ..
                   [=[non-table/non-numeric data.]=])
            if not (k == "THETA_S") then
                assert(
                    not (last_surface[k] == nil),
                    [=[load_last() could not data in specified file for ]=] .. k
                )
                -- Should check that table only contains numeric data
                -- TODO: finish this
                new_initial[k] = htopography_grid.dcopy(last_surface[k])
            end
        end
    end

    new_initial.THETA_S = grid.interp_cubic(last_surface.PHI_S_MID,
                                            last_surface.THETA_S_MID,
                                            _phi_s)

    new_initial.THETA_S[1] = 0 -- given

     new_initial.THETA_S,
    new_initial.BETA_SUB,
       new_initial.N_SUB,
     new_initial.PHI_SUB,
    new_initial.DPHI_SUB = grid.new_grid_from_tau_mid_theta(
                               _phi_s,
                               last_surface.TAU_S_MID,
                               new_initial.THETA_S,
                               n,
                               new_initial.CLUSTER_VAL,
                               new_initial.HOMOTOPY_S,
                               phi_func
                           )

    local _, is_free, is_continued = init_tables_fn(n,
                                                    #(initial_fs.N_SUB))

    DEBUG_OUTPUT(prefix, [=[    load_last()]=], '\n')
    DEBUG_OUTPUT(prefix, [=[        Loaded surface:]=], '\n')
    if VERBOSE then
        surface.display(prefix .. '        ',
                        new_initial,
                        output_data)
    end

    DEBUG_OUTPUT(prefix, [=[    load_last()]=], '\n')
    DEBUG_OUTPUT(prefix, [=[        Creating new surface object.]=], '\n')

    -- TODO: constructor
    local initial = surface.new(new_initial,
                                is_free,
                                is_continued,
                                surface_residual,
                                surface_output,
                                {'A','B'})

    local DIRECTION = continuator.FORWARD

    DEBUG_OUTPUT(prefix, [=[    load_last()]=], '\n')
    DEBUG_OUTPUT(prefix, [=[        Solve via Newton's method.]=], '\n')

    local result = initial.solve_and_reset()

    assert(result == continuator.CONVERGED,
           [=[load_last() could not solve for loaded parameters.]=])
    surf = initial.expand()

    DEBUG_OUTPUT(prefix, [=[    load_last()]=], '\n')
    DEBUG_OUTPUT(prefix, [=[        Converged surface.]=], '\n')
    if VERBOSE then
        surface.display(prefix .. '        ',
                        surf,
                        output_data)
    end

    return initial, is_free, is_continued

end

return {
    load_last = load_last
}

