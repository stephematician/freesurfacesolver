--[[ file: surface.lua

  Some helpers for free surface continuation, which interfaces with continuator
  package.

  Author        : Stephen Wade
  Date created  : ?
  Last modified : 6/2/2013
]]--

require "continuator"

module("surface", package.seeall)

local dcopy
local create_file
local append_file

new = function(in_initial,
               in_is_free,
          in_is_continued, 
              in_residual,
        in_compute_output)

  r_table = {}
  
  r_table.initial = dcopy(in_initial)
  r_table.is_free = dcopy(in_is_free)
  r_table.is_continued = dcopy(in_is_continued)
  
  r_table.residual = in_residual
  
  r_table.compute_output = in_compute_output
  
  r_table.unknowns, r_table.ind_table = pack_vec(in_initial, in_is_free, in_is_continued)
  
  r_table.tangent = {}
  for k, v in pairs(r_table.unknowns) do
    r_table.tangent[k] = 0
  end
  
  t_vec = r_table.residual(r_table.unknowns, {r_table})
  if (#t_vec) ~= (#r_table.unknowns - 1) then
    error('Residual function dimension is ' .. #t_vec .. 
          ' which is != (' .. (#r_table.unknowns) .. ' - 1)')
  end
  
  r_table.append_file = function(in_surface, filename)
    append_file(filename,
                in_surface.unknowns,
                in_surface.initial,
                in_surface.is_free,
                in_surface.is_continued,
                in_surface.ind_table,
                in_surface.compute_output)
  end

  r_table.create_file =
  function(in_surface, filename)
    create_file(filename)
  end
  
  r_table.pcc = continuator.new(in_residual)
  
  r_table.progress = function(in_surface, stepsize, direction, ...)
    in_surface:set_stepsize(stepsize)
    in_surface:set_direction(direction)

    c, t, pr = in_surface.pcc:progress(in_surface.unknowns,
                                        in_surface.tangent,
                                      in_surface.direction,
                                                       arg)
    for k, v in pairs(c) do
      in_surface.unknowns[k] = v
    end
    for k, v in pairs(t) do
      in_surface.tangent[k] = v
    end
    
    return pr
   
  end

  r_table.stepcount = 0
  r_table.stepsize = 0
  r_table.max_stepsize = 0
  r_table.length_vars = {}

  r_table.test_switch = function(in_surface, surf2, surf1, old_length_vars)

    local dp = 0 
    for k = 1, #in_surface.length_vars do
      dp = dp + (surf2[in_surface.length_vars[k]] - surf1[in_surface.length_vars[k]]) *
                ((old_length_vars[k])[2] - (old_length_vars[k])[1])
    end

    return (not (dp > 0))

  end


  r_table.reset_interp = function(in_surface)

    in_surface.stepcount = 0

    in_surface.prev_length_vars = nil

    in_surface.prev_stepsize = nil

    in_surface.prev_finalval = nil

  end

  r_table.set_length_vars = function(in_surface, new_length_vars)

    for k = 1, #new_length_vars do
      assert(in_surface.initial[new_length_vars[k]] ~= nil,
             'Cannot find selected variable(s) for steplength adjustment')
    end

    in_surface.length_vars = {}
    for k = 1, #new_length_vars do
      in_surface.length_vars[k] = new_length_vars[k]
    end

  end

  r_table.set_direction = function(in_surface, new_direction)
   
     assert((new_direction == continuator.BACKWARD) or 
            (new_direction == continuator.FORWARD),
            'Invalid direction')

    in_surface.direction = new_direction

  end

  r_table.get_prev_finalval = function(in_surface)

    return dcopy(in_surface.prev_finalval)

  end

  r_table.get_prev_stepsize = function(in_surface)

    return dcopy(in_surface.prev_stepsize)

  end

  r_table.get_prev_length_vars = function(in_surface)

    return dcopy(in_surface.prev_length_vars)

  end

  r_table.get_length_vars = function(in_surface)

    return dcopy(in_surface.length_vars)

  end

  r_table.get_stepsize = function(in_surface)

    return dcopy(in_surface.stepsize)

  end

  r_table.get_max_stepsize = function(in_surface)

    return dcopy(in_surface.max_stepsize)

  end

  r_table.get_direction = function(in_surface)
   
    return dcopy(in_surface.direction)

  end

  r_table.set_prev_finalval = function(in_surface, new_prev_finalval)

    in_surface.prev_finalval = dcopy(new_prev_finalval)

  end

  r_table.set_stepsize = function(in_surface, new_stepsize)

    assert(new_stepsize >= 0, 'Cannot set stepsize to be negative')

    in_surface.stepsize = new_stepsize
    in_surface.pcc:set_stepsize(new_stepsize)

  end

  r_table.set_max_stepsize = function(in_surface, new_max_stepsize)

    assert(new_max_stepsize > 0, 'Cannot set a maximum stepsize of 0 or less')

    in_surface.max_stepsize = new_max_stepsize

  end

  r_table.set_prev_stepsize = function(in_surface, new_prev_stepsize)

    in_surface.prev_stepsize = dcopy(new_prev_stepsize)

  end

  r_table.set_prev_length_vars = function(in_surface, new_prev_length_vars)

    in_surface.prev_length_vars = dcopy(new_prev_length_vars)

  end

  r_table.progress_approach = function(in_surface, finalval, converge_rate, ...)


    local csurf = unpack_vec(in_surface.unknowns,
                              in_surface.initial,
                              in_surface.is_free,
                         in_surface.is_continued,
                            in_surface.ind_table)

    local stepsize = 0
    local direction = in_surface:get_direction()

    local msg        = '  progress_approach()\n'
    local filler_str = '  -- '
    
    local cfinalval = in_surface.unknowns[#(in_surface.unknowns)]

    stepsize = in_surface.prev_stepsize[2] * converge_rate * 
               (cfinalval - finalval) / 
               (in_surface.prev_finalval[2] - cfinalval)

    if stepsize < 0 then
      stepsize = -stepsize
      
      if direction == continuator.FORWARD then
        direction = continuator.BACKWARD
      else
        direction = continuator.FORWARD
      end

      in_surface:set_direction(direction)
    end

    if stepsize > in_surface.max_stepsize then
      msg      = msg .. filler_str ..'max stepsize used (interpolated stepsize = ' ..
                 tostring(stepsize) .. ').\n'
      stepsize = in_surface.max_stepsize
    end

    msg = msg .. filler_str .. 'stepsize = ' .. tostring(stepsize) .. '.\n'
    local res = in_surface:progress(stepsize, direction, unpack(arg))

    -- update the continuation history
    if in_surface.prev_length_vars then
      for k = 1, #in_surface.length_vars do
         in_surface.prev_length_vars[k] = { (in_surface.prev_length_vars[k])[2], csurf[in_surface.length_vars[k]] }
      end
    else
      in_surface.prev_length_vars = {}
      for k = 1, #in_surface.length_vars do
        in_surface.prev_length_vars[k] = { nil, csurf[in_surface.length_vars[k]] }
      end
    end

    if in_surface.prev_finalval then
      in_surface.prev_finalval = { in_surface.prev_finalval[2],
                                   cfinalval  }
    else
      in_surface.prev_finalval = { nil,
                                   cfinalval }
    end

    if in_surface.prev_stepsize then
      in_surface.prev_stepsize = { in_surface.prev_stepsize[2], stepsize }
    else
      in_surface.prev_stepsize = { nil, stepsize }
    end

    return res, msg
  end

  r_table.progress_interp = function(in_surface, ds, direction, ...)
    
    local csurf = unpack_vec(in_surface.unknowns,
                              in_surface.initial,
                              in_surface.is_free,
                         in_surface.is_continued,
                            in_surface.ind_table)

    -- Calculate a new step size based on continuation history
    local stepsize = 0

    local msg        = '  progress_interp()\n'
    local filler_str = '  -- '

    local cfinalval = in_surface.unknowns[#(in_surface.unknowns)]

    if in_surface.stepcount == 0 then
      msg      = msg .. filler_str .. 'using initial step-size.\n'
      stepsize = in_surface.stepsize
    elseif in_surface.stepcount == 1 then
      -- Use linear 'extrapolation' formula
      msg = msg .. filler_str .. '1st order approximations to derivative available.\n'
      msg = msg .. filler_str .. 'using 1st order truncated T.S. step-size adjustment.\n'
      local ds1 = 0
      local ds2 = ds

      for k = 1, #in_surface.length_vars do
        ds1 = ds1 + math.pow(csurf[in_surface.length_vars[k]] - (in_surface.prev_length_vars[k])[2], 2)
      end
      ds1 = math.sqrt(ds1)

      local step1 = in_surface.prev_stepsize[2]

      stepsize = step1 * ds2 / ds1

    else
      -- Use quadratic 'extrapolation' formula
      msg = msg .. filler_str .. '2nd order approximations to derivatives available.\n'
      local ds1 = 0
      local ds2 = 0
      local ds3 = ds

      for k = 1, #in_surface.length_vars do
        ds1 = ds1 + math.pow((in_surface.prev_length_vars[k])[2] - (in_surface.prev_length_vars[k])[1], 2)
        ds2 = ds2 + math.pow(csurf[in_surface.length_vars[k]] - (in_surface.prev_length_vars[k])[2], 2)
      end
      ds1 = math.sqrt(ds1)
      ds2 = math.sqrt(ds2)

      local step1 = in_surface.prev_stepsize[1]
      local step2 = in_surface.prev_stepsize[2]

      -- This has been changed to use a Taylor series expansion of 's' (arclength)
      -- where the 1st & 2nd derivative of s is approximated at x = x_2 via LIP of
      -- either 1st or 2nd order, depending on whether we get sensible results.
      local sd  = ((ds1+ds2) * (2*step2+step1) / (step2*(step1+step2))) -
                  (ds1 * (step1+step2) / (step1*step2))
      local sdd = (ds1+ds2) / (step2*(step1+step2)) - ds1 / (step1*step2)
      sdd = sdd * 2

      local a = sdd / 2
      local b = sd
      local c = -ds3

      local use_second = false
      if b*b - 4*a*c > 0 then
        msg        = msg .. filler_str .. 'using 2nd order truncated T.S. step-size adjustment.\n'
        use_second = true
        stepsize   = ((-b + math.sqrt(b*b - 4*a*c))/(2*a))
      end

      if (not use_second) or (stepsize < 0) then
        msg      = msg .. filler_str .. 'warning >> resorted to 1st order truncated T.S.\n'
        msg      = msg .. filler_str .. '        >> may require smaller stepsize.\n'
        -- We have also abandoned the second order approximation to the derivative.
        stepsize = step2 * ds3 / ds2
      end

    end

    if stepsize > in_surface.max_stepsize then
      msg      = msg .. filler_str ..'max stepsize used (interpolated stepsize = ' ..
                 tostring(stepsize) .. ').\n'
      stepsize = in_surface.max_stepsize
    end

    msg = msg .. filler_str .. 'stepsize = ' .. tostring(stepsize) .. '.\n'

    in_surface.stepcount = in_surface.stepcount + 1
    local res = in_surface:progress(stepsize, direction, unpack(arg))

    -- update the continuation history
    if in_surface.prev_length_vars then
      for k = 1, #in_surface.length_vars do
         in_surface.prev_length_vars[k] = { (in_surface.prev_length_vars[k])[2], csurf[in_surface.length_vars[k]] }
      end
    else
      in_surface.prev_length_vars = {}
      for k = 1, #in_surface.length_vars do
        in_surface.prev_length_vars[k] = { nil, csurf[in_surface.length_vars[k]] }
      end
    end

    if in_surface.prev_finalval then
      in_surface.prev_finalval = { in_surface.prev_finalval[2], cfinalval }
    else
      in_surface.prev_finalval = { nil, cfinalval }
    end

    if in_surface.prev_stepsize then
      in_surface.prev_stepsize = { in_surface.prev_stepsize[2], stepsize }
    else
      in_surface.prev_stepsize = { nil, stepsize }
    end

    return res, msg

  end

  return r_table
end

-- Taken from lua-users wiki
dcopy = function(object)
  local lookup_table = {}

  local function _copy(object)
    if type(object) ~= "table" then
      return object
    elseif lookup_table[object] then
      return lookup_table[object]
    end
    
    local new_table = {}
    lookup_table[object] = new_table
    
    for index, value in pairs(object) do
      new_table[_copy(index)] = _copy(value)
    end
    
    return setmetatable(new_table, getmetatable(object))
  end

  return _copy(object)
end

pack_vec = function(current, is_free, is_continued)
--[[ This function selects the 'free' parts of a vector
     and arranges them in sequence.
--]]
  packed = {};
  ind_table = {}

  continued_val = nil;
  continued_spec = false;

  i = 1;
  for k, cvec in pairs(current) do
  
    if type(cvec) == "number" then
      
      if is_free[k] then
        if not is_continued[k] then
          packed[i] = cvec
          ind_table[k] = i
          
          i = i + 1;
        else -- continued
          if continued_spec then
            error([=[ Too many variables selected for continuation ]=])
          end
          continued_spec = true
          continued_value = cvec
          ind_table[k] = 0
        end
      else -- not free
        if is_continued[k] then
          error([=[ Non-free variable selected for continuation ]=])
        end
      end
      
    elseif type(cvec) == "table" then
      ind_table[k] = {}
      
      for kk, cnum in pairs(cvec) do
      
        if type(cnum) == "number" then
          if (is_free[k])[kk] then
            if not (is_continued[k])[kk] then
              packed[i] = cnum
              ind_table[k][kk] = i
              
              i = i + 1
            else -- continued
              if continued_spec then
                error([=[ Too many variables selected for continuation ]=])
              end
              continued_spec = true
              continued_value = cnum
              ind_table[k][kk] = 0
            end
          else -- not free
            if (is_continued[k])[kk] then
              error([=[ Non-free variable selected for continuation ]=])
            end
          end
        else -- type cnum not == number
          error([=[ Expected a table containing numbers and tables of numbers
                    for pack vector]=])
        end
        
      end -- Iteration over the table current[cvec]
      
    else -- type of cvec not == number or table 
      error([=[ Expected a table containing numbers and tables of numbers for
             pack vector]=])
    end
    
  end -- Iteration over pairs(current)
  
  packed[i] = continued_value

  for k, v in pairs(ind_table) do
    if type(v) == "table" then
      for kk, vv in pairs(v) do
        if vv == 0 then
          ind_table[k][kk] = i
        end
      end
    else
      if v == 0 then
        ind_table[k] = i
      end
    end
  end

  return packed, ind_table
  
end

unpack_vec = function(packed, initial, is_free, is_continued, ind_table)
--[[ Comment?
--]]
  current = {}

  c_specified = false

  for k, cvec in pairs(initial) do
    
    if type(cvec) == "number" then
      if is_free[k] then
        if not is_continued[k] then
          i = ind_table[k]
          current[k] = packed[i]
        else -- continued
          if c_specified then
            error([=[ Too many variables selected for continuation ]=])
          end
          i = ind_table[k]
          current[k] = packed[i]
          --current[k] = packed[#packed]
          c_specified = true
        end
      else -- not free
        if is_continued[k] then
          error([=[ Non-free variable selected for continuation ]=]) 
        end
        current[k] = cvec;
      end
  
    elseif type(cvec) == "table" then
      current[k] = {}
      for kk, cnum in pairs(cvec) do
        if (is_free[k])[kk] then
          if not (is_continued[k])[kk] then 
            i = (ind_table[k])[kk];
            (current[k])[kk] = packed[i]
          else -- continued
            if c_specified then
              error([=[ Too many variables selected for continuation ]=])
            end
            i = (ind_table[k])[kk];
             (current[k])[kk] = packed[i]
           -- (current[k])[kk] = packed[#packed]
            c_specified = true   
          end
        else -- not free
          if (is_continued[k])[kk] then
            error([=[ Non-free variable selected for continuation ]=]) 
          end
          (current[k])[kk] = cnum;
        end
        
      end -- Iteration over the table initial[cvec]
      
    else -- type of cvec not == number or table
      error([=[ Expected a table containing numbers and tables of numbers for
             pack vector]=]);
    end
    
  end -- Iteration over pairs(initial)
  
  return current
  
end

--create_file = function(filename, delta)
create_file = function(filename)

  f, l = string.find(filename, ".", 1, true)

  if l then
    ff, ll = string.find(filename, "/", 1, true)
    if ll then --?
      fn_name = string.sub(filename, ll+1, l-1)
    else
      fn_name = string.sub(filename, 1, l-1)
    end
  else
    ff, ll = string.find(filename, "/", 1, true)
    if ll then
      fn_name = string.sub(filename, ll+1, string.length(filename))
      else
    fn_name = string.sub(filename, 1, string.length(filename))
    end
  end

  of = io.open(filename, "w")
  of:write(string.format("%% %08d", 0) .. "\n\n")
  of:write("function [ o ] = " .. fn_name .. "()\n\n")
  of:write("\nend\n")
  of:close()
    
end

append_file = function(filename,
                         packed,
                        initial,
                        is_free,
                   is_continued,
                      ind_table,
                 compute_output)--,

  of = io.open(filename, "r+")

  f, l = string.find(filename, ".", 1, true)

  if l then
    ff, ll = string.find(filename, "/", 1, true)
    if ll then
      fn_name = string.sub(filename, ll+1, l-1)
    else
      fn_name = string.sub(filename, 1, l-1)
    end
  else
    ff, ll = string.find(filename, "/", 1, true)
    if ll then
      fn_name = string.sub(filename, ll+1, string.length(filename))
      else
    fn_name = string.sub(filename, 1, string.length(filename))
    end
  end
  
  of:seek("set", string.len("%% "))
  n = of:read("*n")
  n = n + 1
  of:seek("set", 0)
  of:write(string.format("%% %08d", n) .. "\n\n")
  of:write("function [ o ] = " .. fn_name .. "()\n\n")

  of:close()
                                
  temp = unpack_vec(packed, initial, is_free, is_continued, ind_table)
  
  out_table = compute_output(temp)
  
  of = io.open(filename, "r+")
  of:seek("end", -string.len("end\n"))
    
  for k, v in pairs(out_table) do
    if type(v) == "number" then
    
      of:write("  o." .. k .. "{" .. n .. "} = " .. v .. ";\n")
      
    elseif type(v) == "table" then
    
      of:write("  o." .. k .. "{" .. n .. "} = [")
      
      for kk, vv in pairs(v) do

        of:write(vv .. " ")
      end
      
      of:write("];\n")
      
    else -- wtf
      error("Output computed of surface should be vectors and scalars only")
    end
  end

  of:write("\nend\n")
  
  of:close()
  
end
