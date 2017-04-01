--[[
Compute continuation curve of box shaped disturbance with fixed width, and fixed
Froude number.

Performs two or three continuations. The first is from uniform stream, and it
continues clockwise around until the velocity at the crest gets too small, or
it hits double hump solutions with a small box-height. It then performs a
continuation from the original uniform stream to the 'dip' like solutions.
If no double hump solutions were found in the first continuation, then this
finds a double hump solution and then performs a clockwise continuation.

rev03 - homotopy available, symmetric formulation, with dbeta dependent
        upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 7/12/2012
Modified - 25/3/2013

History:
25/3/2013 Updated documentation.

8/3/2013  Added in 'emailme', developed ability to deal with the branches being
          either connected or not connected.

4/3/2013  Modified to deal with new script layout in rev67 of code

25/2/2013 Computation for u = 0.2 at the crest, removed box_height > 0 stuff.

16/2/2013 Converted to new way of calling continuation method progress_interp().

7/12/2012 Created original script.
]]--

require "surface"
require "pressure_residuals"
require "pressure_grid"

DO_US1 = true
DO_US2 = true
DO_DBL = true

EMAILME = true
OVERWRITE = true

-- Grid
N_HI_RES = 600
N_LO_RES = 200

CLUSTER_DPHI_VAL = 0.01

-- Domain and topography
B       = 2.80
FROUDE  = 1.32
PHI_LIM = 10

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

FSTR = string.format('%3u',round(FROUDE*100))
BSTR = string.format('%3u',round(B*100))

-- Starting and finishing
UNIFORM_A    = -0.05
STOKES_U_MIN = 0.05
LIMIT_DBL_A  = -0.0005*B
LARGEST_AREA = 8 -- 10 or 15?

in_script_us  = 'input/pressure_uniform_init.lua'
in_script_dbl = 'input/pressure_double_init.lua'

out_file_us1 = 'output/pressure_F' .. FSTR .. '_B' .. BSTR .. '.m'
out_file_us2 = 'output/pressure_F' .. FSTR .. '_DIP_B' .. BSTR .. '.m'
out_file_dbl = 'output/pressure_F' .. FSTR .. '_DOUBLE_B' .. BSTR .. '.m'

found_double = false

if DO_US1 or DO_US2 then
  dofile(in_script_us)
end

-- Email notification sent when computation is completed.
if EMAILME then
  require "luarocks.loader"
  smtp = require("socket.smtp")

  from = "<stephen.wade@adelaide.edu.au>"
  smtpserver = "smtp.adelaide.edu.au"

  rcpt = {
    "<stephematician@gmail.com>"
  }
end

mesgt = {
  headers = {
    to = "Stephen Wade <stephematician@gmail.com>",
    from = "Stephen Wade <stephen.wade@adelaide.edu.au>",
    subject = "pressure calculation" 
  },
  body = "Body:\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}


if DO_US1 then

freesurface = surface.new(    hi_res_unif,
                         free_hi_res_unif,
                    continued_hi_res_unif,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)


--[[ Output file generation ]]--

of = io.open(out_file_us1, "r")

if not of then
   freesurface:create_file(out_file_us1)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_us1)
    freesurface:create_file(out_file_us1) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

num_continuate = 400
default_ds = 0.01
ds = default_ds

stepsize = 0.01
max_stepsize = 2.5

freesurface:set_max_stepsize(max_stepsize)
freesurface:set_stepsize(stepsize)

direction = continuator.BACKWARD

io.write('max iterations : ' .. tostring(num_continuate) .. '\n')
io.write('     step size : ' .. tostring(stepsize) ..'\n')
io.flush()

initial_step = true
res = nil

for i = 1, num_continuate, 1 do
  local outstr = 'Iteration : ' .. tostring(i-1)
  local filstr = '-'
  io.write(outstr .. '\n')
  io.write(filstr:rep(outstr:len()) .. '\n')

  local csurf = surface.unpack_vec(freesurface.unknowns,
                                    freesurface.initial,
                                    freesurface.is_free,
                               freesurface.is_continued,
                                  freesurface.ind_table)

  io.write('  Parameter summary :\n')  
  io.write('        A : ', csurf.A, '\n')
  io.write('     ETA0 : ', csurf.ETA0, '\n')
  io.write('    U_MIN : ', csurf.U_MIN, '\n')
  io.flush()

  if res == continuator.BIFURCATION then
    io.write('  Bifurcation detected.\n')
    io.flush()
    if direction == continuator.FORWARD then
      direction = continuator.BACKWARD
    else
      direction = continuator.FORWARD
    end
    
  end

  if res == continuator.DIVERGED then
    mesgt.body = mesgt.body .. 
                 'Finished due to lack of convergence in continuation method.\n'
    io.write('  Divergence in continuation.\n')
    io.flush()
    break
  end

  if res == continuator.ERROR then
    mesgt.body = mesgt.body .. 
                 'Finished due to error in continuation method.\n'
    error('Error detected in continuation.')
    break
  end

  freesurface:append_file(out_file_us1)

  if csurf.U_MIN < STOKES_U_MIN then
    mesgt.body = mesgt.body .. 
                 'Finished first continuation branch at a Stokes configuration.\n' ..
                 'Took ' .. i .. ' iterations.\n'
    break
  end

  local v = csurf.V_S
  local is_double = false
  for j = 2, #v do
     if (v[j] * v[j-1] < 0) and (v[j-1] > 0) then
      is_double = true
      break
    end
  end

  if is_double and csurf.A > LIMIT_DBL_A then
    found_double = true
    mesgt.body = mesgt.body .. 
                 'Finished first continuation branch at a double hump solution.\n' ..
                 'Took ' .. i .. ' iterations.\n'
    break
  end

  if csurf.U_MIN < 0.15 then
    ds = default_ds*0.5
  elseif csurf.U_MIN < 0.40 then
    ds = default_ds*0.8
  else
    ds = default_ds
  end

  local progress_start = os.clock()
  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  local progress_finish = os.clock()
  io.write(msg)
  io.write('    Continuation time = ' .. tostring(progress_finish - progress_start) ..'\n')
  io.flush()

  if initial_step then
    local nsurf = surface.unpack_vec(freesurface.unknowns,
                                      freesurface.initial,
                                      freesurface.is_free,
                                 freesurface.is_continued,
                                    freesurface.ind_table)

    if nsurf.A < csurf.A then
      freesurface:reset_interp()
      io.write('  Started in wrong direction, attempting to switch.\n')
      if direction == continuator.FORWARD then
        direction = continuator.BACKWARD
      else
        direction = continuator.FORWARD
      end
    else
      initial_step = false
    end
  else
    local regrid_start = os.clock()
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage('collect')

end

end -- DO_US1

if DO_US2 then

--[[
Find the dip-like solutions? Then continue on to a pseudo-limiting configuration
with large aspect ratio of disturbance?
--]]

freesurface = surface.new(    hi_res_unif,
                         free_hi_res_unif,
                    continued_hi_res_unif,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

of = io.open(out_file_us2, "r")

if not of then
   freesurface:create_file(out_file_us2)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_us2)
    freesurface:create_file(out_file_us2) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

default_ds = 0.015
num_continuate = 400
ds = default_ds

stepsize = 0.01
max_stepsize = 2.5

freesurface:set_max_stepsize(max_stepsize)
freesurface:set_stepsize(stepsize)

direction = continuator.BACKWARD

io.write('max iterations : ' .. tostring(num_continuate) .. '\n')
io.write('     step size : ' .. tostring(stepsize) ..'\n')
io.flush()

initial_step = true
res = nil

for i = 1, num_continuate, 1 do
  local outstr = 'Iteration : ' .. tostring(i-1)
  local filstr = '-'
  io.write(outstr .. '\n')
  io.write(filstr:rep(outstr:len()) .. '\n')

  local csurf = surface.unpack_vec(freesurface.unknowns,
                                    freesurface.initial,
                                    freesurface.is_free,
                               freesurface.is_continued,
                                  freesurface.ind_table)

  io.write('  Parameter summary :\n')  
  io.write('        A : ', csurf.A, '\n')
  io.write('     ETA0 : ', csurf.ETA0, '\n')
  io.write('    U_MIN : ', csurf.U_MIN, '\n')
  io.flush()

  if res == continuator.BIFURCATION then
    io.write('  Bifurcation detected.\n')
    io.flush()
    if direction == continuator.FORWARD then
      direction = continuator.BACKWARD
    else
      direction = continuator.FORWARD
    end
    
  end

  if res == continuator.DIVERGED then
    mesgt.body = mesgt.body .. 
                 'Finished due to lack of convergence in continuation method.\n'
    io.write('  Divergence in continuation.\n')
    io.flush()
    break
  end

  if res == continuator.ERROR then
    mesgt.body = mesgt.body .. 
                 'Finished due to error in continuation method.\n'
    error('Error detected in continuation.')
    break
  end

  freesurface:append_file(out_file_us2)

  if math.abs(csurf.A * csurf.B) > LARGEST_AREA then
    mesgt.body = mesgt.body .. 
                 'Finished dip-like solutions with reasonably large dip.\n' ..
                 'Took ' .. i .. ' iterations.\n'
    break
  end

  local progress_start = os.clock()
  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  local progress_finish = os.clock()
  io.write(msg)
  io.write('    Continuation time = ' .. tostring(progress_finish - progress_start) ..'\n')
  io.flush()

  if initial_step then
    local nsurf = surface.unpack_vec(freesurface.unknowns,
                                      freesurface.initial,
                                      freesurface.is_free,
                                 freesurface.is_continued,
                                    freesurface.ind_table)

    if nsurf.A > csurf.A then
      freesurface:reset_interp()
      io.write('  Started in wrong direction, attempting to switch.\n')
      if direction == continuator.FORWARD then
        direction = continuator.BACKWARD
      else
        direction = continuator.FORWARD
      end
    else
      initial_step = false
    end
  else
    local regrid_start = os.clock()
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage("collect")

end

end -- if DO_US2

if DO_DBL then

if not found_double then
dofile(in_script_dbl)

freesurface = surface.new(  hi_res_sol,
                       free_hi_res_sol,
                  continued_hi_res_sol,
        pressure_residuals.residual_hs,
  pressure_residuals.compute_output_hs)

--[[ Output file generation --]]

of = io.open(out_file_dbl, "r")

if not of then
   freesurface:create_file(out_file_dbl)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_dbl)
    freesurface:create_file(out_file_dbl) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

num_continuate = 400
default_ds = 0.015
ds = default_ds

stepsize = 0.001
max_stepsize = 2.5

freesurface:set_max_stepsize(max_stepsize)
freesurface:set_stepsize(stepsize)

direction = continuator.BACKWARD

io.write('max iterations : ' .. tostring(num_continuate) .. '\n')
io.write('     step size : ' .. tostring(stepsize) ..'\n')
io.flush()

initial_step = true
res = nil

for i = 1, num_continuate, 1 do
  local outstr = 'Iteration : ' .. tostring(i-1)
  local filstr = '-'
  io.write(outstr .. '\n')
  io.write(filstr:rep(outstr:len()) .. '\n')

  local csurf = surface.unpack_vec(freesurface.unknowns,
                                    freesurface.initial,
                                    freesurface.is_free,
                               freesurface.is_continued,
                                  freesurface.ind_table)
  
  io.write('  Parameter summary :\n')  
  io.write('        A : ', csurf.A, '\n')
  io.write('     ETA0 : ', csurf.ETA0, '\n')
  io.write('    U_MIN : ', csurf.U_MIN, '\n')
  io.flush()

  if res == continuator.BIFURCATION then
    io.write('  Bifurcation detected.\n')
    io.flush()
    if direction == continuator.FORWARD then
      direction = continuator.BACKWARD
    else
      direction = continuator.FORWARD
    end
    
  end

  if res == continuator.DIVERGED then
    mesgt.body = mesgt.body .. 
                 'Finished due to lack of convergence in continuation method.\n'
    io.write('  Divergence in continuation.\n')
    io.flush()
    break
  end

  if res == continuator.ERROR then
    mesgt.body = mesgt.body .. 
                 'Finished due to error in continuation method.\n'
    error('Error detected in continuation.')
    break
  end

  freesurface:append_file(out_file_dbl)

  if math.abs(csurf.A * csurf.B) > LARGEST_AREA then
    mesgt.body = mesgt.body .. 
                 'Finished doube-hump solutions with reasonably large trench.\n' ..
                 'Took ' .. i .. ' iterations.\n'
    break
  end

  if csurf.U_MIN < 0.15 then
    ds = default_ds * 0.25 -- 0.0025
  elseif csurf.U_MIN < 0.40 then
    ds = default_ds * 0.5 -- 0.0050
  else
    ds = default_ds
  end

  local progress_start = os.clock()
  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  local progress_finish = os.clock()
  io.write(msg)
  io.write('    Continuation time = ' .. tostring(progress_finish - progress_start) ..'\n')
  io.flush()

  if initial_step then
    local nsurf = surface.unpack_vec(freesurface.unknowns,
                                      freesurface.initial,
                                      freesurface.is_free,
                                 freesurface.is_continued,
                                    freesurface.ind_table)

    if nsurf.U_MIN < csurf.U_MIN then
      freesurface:reset_interp()
      io.write('  Started in wrong direction, attempting to switch.\n')
      if direction == continuator.FORWARD then
        direction = continuator.BACKWARD
      else
        direction = continuator.FORWARD
      end
    else
      initial_step = false
    end
  else
    local regrid_start = os.clock()
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage("collect")

end

end -- if not found_double

end -- if DO_DBL

io.write('Finished continuation for all solutions.\n')
io.flush()

if EMAILME then
  r, e = smtp.send{
    from = from,
    rcpt = rcpt, 
    source = smtp.message(mesgt),
    server = smtpserver
  }
else
  print(mesgt.headers.subject)
  print(mesgt.body)
end
