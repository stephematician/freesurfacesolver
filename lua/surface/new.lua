local unknowns    = require "surface.unknowns"
local iosurf      = require "surface.io"
local expansions  = require "surface.expansions"
local continuator = require "continuator"

--[[
--  Create a new surface object for continuation
--
--    - _initial Table of initial values (numeric or vector) of surface
--    - _is_free Table of boolean indicating free or fixed for each variable
--    - _is_continued Table with one boolean value indicating the continued
--                    variable
--    - _residual Function which takes two arguments; the unknowns and a
--                function which expands the unknowns into a table like
--                _initial. Return value is a vector of residuals of the system
--                of equations
--    - _compute_output Function which takes one argument; a table in the same
--                      format as _initial which has the current values of the
--                      surface
--    - _curve_domain A table which contains keys corresponding to numeric
--                    values in _initial which will be used to determine step
--                    sizes so that the surface.progress_interp() function
--                    produces surfaces which are evenly spaced in the domain
--                    specified (by the keys)
--    - _max_step The maximum initial step size for continuation
-- TODO: update this documentation
--]]

local new

new = function(_initial,
               _is_free,
               _is_continued,
               _residual,
               _compute_output,
               _curve_domain,
               _first_step,
               -- Optional arguments
               __max_step,
               __direction,
               __step_history,
               __curve_history,
               __approach_history)


    --[[--------------------------------------------------------------
    --                                         Forward declarations --
    ----------------------------------------------------------------]]
    
    --[[
    --  Copy constructor.
    --  
    --  If no arguments are supplied (default behaviour) then the whole surface
    --  object is copied.
    --
    --  If a surface description _is_ supplied, then a new surface is created
    --  with the same residual and compute output functions, but none of the
    --  continuation history etc.
    --]]
    local copy
    --[[
    -- Expand the unknowns (a vector/array) into a table with keys.
    --]]
    local expand
    --[[
    -- Create a file for output.
    --]]
    local create_file
    --[[
    -- Append the table produced by the compute_output() member function to a
    -- matlab .m file.
    --]]
    local append_file
    --[[
    -- Accessors for step data.
    --]]
    local get_step, set_step
    --[[
    -- Accessors for direction data.
    --]]
    local get_direction, set_direction
    --[[
    -- Accessor for maximum step length parameter.
    --]]
    local get_max_step, set_max_step
    --[[
    -- Accessor for the continuation curve 'domain' variables, i.e. the
    -- domain in which we'd like the continuation curve to be presented.
    -- Typically this is the amplitude of forcing-eta0 plane.
    --]]
    local get_curve_domain, set_curve_domain
    --[[
    -- Accessor for step (length) history, the continuation curve history, and
    -- the history of the approach to a specific variable.
    --]]
    local get_step_history, get_curve_history, get_approach_history
    --[[
    -- Flush the step, curve and approach histories.
    --]]
    local flush_history
    --[[
    -- Accessor for surface description.
    --]]
    local get_is_free, get_is_continued
    --[[
    -- Perform one step of the continuation algorithm.
    --]]
    local progress 
    --[[
    -- Test if the curve traced from two given surfaces is in the same direction
    -- as the curve traced in the last continuation step.
    --]]
    local test_switch
    --[[
    -- Attempt to solve for the current value of the continuation parameter
    --]]
    local solve_and_reset
    --[[
    -- Attempt to approach a value of the continuation parameter
    --
    --   - approach_value desired continuation parameter value
    --   - converge_rate the coefficient of the step size used to control
    --                   the approach
    --]]
    local progress_approach
    --[[
    -- Attempt to approach a value of the continuation parameter
    --
    --   - ds desired arclength in the plotting domain of step
    --   - _direction direction to perform continuation
    --]]
    local progress_interp

    --[[--------------------------------------------------------------
    --                                      Default argument values --
    ----------------------------------------------------------------]] 
    local _direction, _max_step = __direction, __max_step
    
    if not _direction then
        _direction = continuator.FORWARD
    end
    if not _max_step then
        _max_step = 1/0
    end
    
    local _step_history     = __step_history
    local _curve_history    = __curve_history
    local _approach_history = __approach_history

    if not _step_history then
        _step_history = stack.new(2)
    end
    if not _curve_history then
        _curve_history = stack.new(3)
    end
    if not _approach_history then
        _approach_history = stack.new(3)
    end

    --[[--------------------------------------------------------------
    --                                  Private data initialisation --
    ----------------------------------------------------------------]]
    local self = {
        -- Surface description - should have N unknowns
        initial     = iosurf.deepcopy(_initial),
        is_free      = iosurf.deepcopy(_is_free),
        is_continued = iosurf.deepcopy(_is_continued),
        -- Packed surface description
        unknowns, index_map = unknowns.pack(_initial,
                                            _is_free,
                                            _is_continued),
        -- Residual for system of N-1 equations
        residual       = _residual,
        compute_output = _compute_output,
        -- Tangent vector in N-dimensionsal 'unknowns' space
        tangent = {},
        -- Continuation plotting/curve parameters
        curve_domain = set_curve_domain(_curve_domain),
        max_step     = set_max_step(_max_step),
        -- Continuation plotting/curve data
        step_history     = _step_history,
        curve_history    = _curve_history,
        approach_history = _approach_history
    }

    if not __curve_history then
        local _new_curve = {}
        for k = 1, #self.curve_domain do
            _new_curve[k] = self.initial[self.curve_domain[k]]
        end
        self.curve_history.push(_new_curve)
    end

    if not __approach_history then
        self.approach_history.push(self.unknowns[#self.unknowns])
    end

    for k = 1, #self.unknowns do
        self.tangent[k] = 0
    end

    local _v = _residual(self.unknowns,
                         expand)
    if not ((#_v) == (#self.unknowns - 1)) then
        error([=[Dimension of residual function (]=] ..  #_v .. 
              [=[) should be ]=] ..  #self.unknowns - 1 ..  [=[.]=])
    end

    -- initialise the continuator object
    self.pcc = continuator.new(_residual)
    set_step(_first_step)

    --[[--------------------------------------------------------------
    --                                                     Closures --
    ----------------------------------------------------------------]]

 
    local copy = function(_initial,
                          _is_free,
                          _is_continued)
        if _initial then
            return new(_initial,
                       _is_free,
                       _is_continued,
                       self.residual,
                       self.compute_output,
                       get_curve_domain(),
                       get_step(),
                       get_max_step())
        else 
            local _surface = new(self.initial,
                                 self.is_free,
                                 self.is_continued,
                                 self.residual,
                                 self.compute_output,
                                 get_curve_domain(),
                                 get_step(),
                                 get_max_step(),
                                 get_direction(),
                                 get_step_history(),
                                 get_curve_history(),
                                 get_approach_history())
            return _surface
        end
    end


    local expand = function(_unknowns)

        if not (_unknowns) then
            _unknowns = self.unknowns
        end
            
        return unknowns.expand(_unknowns,
                               self.initial
                               self.is_free,
                               self.ind_table)

    end


    local create_file = function(filename)
        iosurf.create_file(filename)
    end


    local append_file = function(filename)

        return iosurf.append_file(
                   filename,
                   self.compute_output(
                       unknowns.expand(self.unknowns,
                                       self.initial,
                                       self.is_free,
                                       self.is_continued,
                                       self.ind_table)
                   )
               )

    end


    local get_step = function()
        return self.step
    end


    local set_step = function(_step)
        assert(_step >= 0,
               [=[Step size must be greater than or equal to zero]=])
        self.step = _step
        self.pcc:set_step(_step)
    end


    local get_direction = function()
        return self.direction
    end


    local set_direction = function(_direction)
        assert(_direction == continuator.FORWARD or
                   _direction == continuator.BACKWARD,
               [=[Invalid direction]=])
        self.direction = _direction
    end


    local get_curve_domain = function()
        return iosurf.deepcopy(self.curve_domain)
    end


    local set_curve_domain = function(_curve_domain)
        for k = 1, #_curve_domain do
            assert(self.initial[_curve_domain[k]],
                [=[Selected variable for continuation curve domain not ]=] ..
                [=[found in surface description tables.]=]
            )
            assert(
                type(self.initial[_curve_domain[k]]) == 'number',
                [=[Selected variable for continuation curve domain must ]=] ..
                [=[be a name of a numeric in the surface description.]=]
            )
        end

        self.curve_domain = {}
        for k = 1, #_curve_domain do
            self.curve_domain[k] = _curve_domain[k]
        end
    end

    
    local get_max_step = function()
        return self.max_step
    end
    
    
    local set_max_step = function(_max_step)
        assert(_max_step > 0,
               [=[set_max_step() cannot set maximum step of zero or less.]=])
        self.max_step = _max_step
    end


    local get_step_history = function()
        return iosurf.deepcopy(self.step_history)
    end


    local get_approach_history = function()
        return iosurf.deepcopy(self.approach_history)
    end


    local get_curve_history = function()
        return iosurf.deepcopy(self.curve_history)
    end


    local flush_history = function()

        self.curve_history.flush()
        self.step_history.flush()
        self.approach_history.flush()
        
        -- Initialise the curve and approach histories
        local _new_curve = {}
        for k = 1, #self.curve_domain do
            _new_curve[k] = self.initial[self.curve_domain[k]]
        end

        self.curve_history.push(_new_curve)
        self.approach_history.push(self.unknowns[#self.unknowns])

    end


    local get_is_free = function()

        return iosurf.deepcopy(self.is_free)

    end


    local get_is_continued = function()

        return iosurf.deepcopy(self.is_continued)

    end




    local progress = function(_step,
                              _direction)

        set_step(_step) 
        set_direction(_direction)

        c, t, pr = self.pcc:progress(self.unknowns,
                                     self.tangent,
                                     self.direction,
                                     expand)
        -- Deep copy values
        for k, v in ipairs(c) do
            self.unknowns[k] = v
        end
        for k, v in ipairs(t) do
            self.tangent[k] = v
        end

        return pr

    end


    local test_switch = function(s_a, s_b)

        local dp = 0
        local c_b = self.curve_history.peek(1)
        local c_a = self.curve_history.peek(2)

        for k = 1, #self.curve_domain do
            dp = dp + (s_b[self.curve_domain[k]] - 
                           s_a[self.curve_domain[k]]) *
                       (c_b[k] - c_a[k])
        end

        return (not (dp > 0))

    end


    local solve_and_reset = function()

        c, pr = self.pcc:newton(self.unknowns,
                                expand)
        -- deep copy values
        for k, v in ipairs(c) do
            self.unknowns[k] = v
        end

        flush_history()

        return pr

    end


    local progress_approach = function(approach_value, converge_rate)

        local msg        = '    progress_approach()'
        local filler_str = '\n        - '
        local step = 0

        if (self.approach_history.count >= 2) then
            step = self.step_history.peek(1) * converge_rate * 
                        (approach_value - self.approach_history.peek(1)) / 
                        (self.approach_history.peek(1) - 
                             self.approach_history.peek(2))
        else
            -- use initial step size
            msg = msg .. filler_str ..
                  [=[Previous step sizes not available, using surface ]=] ..
                  [=[value.]=]
            step = get_step()
        end

        local direction = get_direction()
        
        if step < 0 then
            step = -step
      
            if direction == continuator.FORWARD then
                direction = continuator.BACKWARD
            else
                direction = continuator.FORWARD
            end

            msg = msg .. filler_str ..
                  [=[Step size estimated as negative, reversing ]=] ..
                  [=[continuation direction.]=]
 
            set_direction(direction)
        end

        if step > get_max_step() then
            msg = msg .. filler_str ..
                  [=[Maximum step size used instead of interpolated ]=] ..
                  [=[stepsize = ]=] .. tostring(step)
            step = get_max_step()
        else
            msg = msg .. filler_str ..
                  [=[Interpolated step size = ]=] .. tostring(step)
        end

        local res = progress(step, direction)

        -- Update the continuation history with new values
        _surface = expand()

        local new_curve = {}
        for k = 1, #curve_domain do
            new_curve[k] = _surface[self.curve_domain[k]]
        end

        self.curve_history.push(new_curve)

        self.approach_history.push(self.unknowns[#self.unknowns])

        self.step_history.push(step)

        return res, msg

    end

    local progress_interp = function(ds, _direction)
    
        local _surface = expand()

        -- Calculate a new step size based on continuation history
        local step = 0

        local msg        = '    progress_interp()'
        local filler_str = '\n        - '

        if self.step_history.get_count() == 0 then
            msg  = msg .. filler_str .. [=[Using initial step size.]=]
            step = get_step()
        elseif self.approach_history.get_count() == 1 then
            -- Use linear 'extrapolation' formula
            msg = msg .. filler_str .. [=[Using first order truncated T.S.]=] ..
                  [=[for arclength to estimate step size]=]
            step = expansions.first_order_stepsize(self.step_history,
                                                   self.curve_history,
                                                   ds)
        else
            step = expansions.second_order_stepsize(self.step_history,
                                                    self.curve_history,
                                                    ds)
            if step < 0 then
                -- Couldn't estimate step size using second order TS expansion
                step = expansions.first_order_stepsize(self.step_history,
                                                       self.curve_history,
                                                       ds)
                msg = msg .. filler_str .. [=[Using first order truncated ]=] ..
                      [=[T.S. of arclength to estimate step size.]=]

            else
                msg = msg .. filler_str .. [=[Using second order ]=] ..
                      [=[truncated T.S. of arclength to estimate step size.]=]

            end
        end

        if step > get_max_step() then
            msg = msg .. filler_str ..
                  [=[Maximum step size ]=] ..
                  tostring(get_max_step) ..
                  [=[ used instead of interpolated step size = ]=] ..
                  tostring(step) 
            step = get_max_step()
        else
            msg = msg .. filler_str ..
                  [=[Interpolated step size = ]=] .. tostring(step)
        end

        local res = progress(step, _direction)

        -- Update the continuation history with new values
        _surface = expand()

        local new_curve = {}
        for k = 1, #curve_domain do
            new_curve[k] = _surface[self.curve_domain[k]]
        end

        self.curve_history.push(new_curve)

        self.approach_history.push(self.unknowns[#self.unknowns])

        self.step_history.push(step)

        return res, msg

    end
    

    return {
        copy = copy,
        expand = expand,
        create_file = create_file,
        append_file = append_file,
        get_step = get_step,
        set_step = set_step,
        get_direction = get_direction,
        set_direction = set_direction,
        get_max_step = get_max_step,
        set_max_step = set_max_step,
        get_curve_domain = get_curve_domain,
        set_curve_domain = set_curve_domain,
        get_step_history     = get_step_history,
        get_curve_history    = get_curve_history,
        get_approach_history = get_approach_history,
        flush_history        = flush_history,
        get_is_free      = get_is_free,
        get_is_continued = get_is_continued,
        test_switch       = test_switch,
        solve_and_reset   = solve_and_reset,
        progress_approach = progress_approach,
        progress_interp   = progress_interp
    }
   
end

return new

