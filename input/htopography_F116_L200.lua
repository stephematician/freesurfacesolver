--[[
Compute continuation curve of smooth box shaped disturbance with fixed width,
and fixed Froude number.

Box is defined via the use of an approximation to a step function given by
the hyperbolic tangent.

Performs two or three continuations. The first is from uniform stream, and it
continues anti-clockwise around until the velocity at the crest gets too small,
or it hits double hump solutions with a small box-height. It then performs a
continuation from the original uniform stream to the 'dip' like solutions.
If no double hump solutions were found in the first continuation, then this
finds a double hump solution with small box height and then performs a clockwise
continuation.

homotopy available, symmetric formulation, with dbeta dependent
upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 07/12/2012
Modified - 11/08/2014

History:
11/08/2014 Added ability to 'restart' calculations.
29/04/2014 In 'beta' phase, appears to function correctly for F=1.1, needs to
           be tried for F=1.11 - 1.30.
09/01/2014 Adapted from gtopography_F110_B200.lua and new initialisation format
           from the pressure_U* code.

25/03/2013 Updated documentation.

08/03/2013 Added in 'emailme', added code to deal with the branches being
           either connected or not connected.

04/03/2013 Modified to deal with new script layout in rev67 of code

25/02/2013 Computation for u = 0.2 at the crest, removed box_height > 0 stuff.

16/02/2013 Converted to new way of calling continuation method progress_interp().

07/12/2012 Created original script.
--]]

require "surface"
require "htopography_residuals"
require "htopography_grid"
require "htopography_load"

dofile('input/htopography_find_pert_to_unif.lua')
dofile('input/htopography_find_doublelim_from_unif.lua')

DO_US1 = true
DO_US2 = true 
DO_DBL = true
DO_DS2 = true

EMAILME = false
OVERWRITE = false

-- Grid
N_HI_RES = 440
N_LO_RES = 240

CLUSTER_DPHI_VAL = 0.01

-- Domain and topography
L       = 2.00
B       = 10.00 -- ?
FROUDE  = 1.16
PHI_LIM = 10

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

FSTR = string.format('%3u',round(FROUDE*100))
LSTR = string.format('%3u',round(L*100))

-- Starting and finishing
UNIFORM_A    = -0.04
LIMIT_DBL_A  = -0.0005*L
STOKES_U_MIN = 0.08
LARGEST_AREA = 3.0
SMALLEST_AREA = 0.005

in_script_us   = 'input/htopography_uniform_init.lua'
in_script_dbl  = 'input/htopography_double_init.lua'

out_file_us1 = 'output/htopography_F' .. FSTR .. '_US1_L' .. LSTR .. '.m'
out_file_us2 = 'output/htopography_F' .. FSTR .. '_US2_L' .. LSTR .. '.m'
out_file_dbl = 'output/htopography_F' .. FSTR .. '_DBL_L' .. LSTR .. '.m'
out_file_ds2 = 'output/htopography_F' .. FSTR .. '_DS2_L' .. LSTR .. '.m'

found_double = false

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
    subject = "htopography calculation, F = " .. FSTR .. ', L = ' .. LSTR 
  },
  body = "Body:\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}


if DO_US1 then

of = io.open(out_file_us1, "r")

if not of then
  io.write('Creating : ' .. out_file_us1 .. '.\n')
  io.flush()
  dofile(in_script_us)
  IN_SCRIPT_US = true
  freesurface:create_file(out_file_us1)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_us1 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_us1)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_us1 .. '.\n')
    io.flush()
    dofile(in_script_us)
    freesurface:create_file(out_file_us1)
    IN_SCRIPT_US = true
  else
    io.write('Loading : ' .. out_file_us1 .. '.\n')
    io.flush()
              hi_res_unif,
         free_hi_res_unif, 
    continued_hi_res_unif = htopography_load.load_last(out_file_us1,
                                                pert_to_unif_tables,
                                                           N_HI_RES)
    IN_SCRIPT_US = false
  end
end

-- Experimenting here.
continued_hi_res_unif.A = false
continued_hi_res_unif.D = true

freesurface = surface.new(    hi_res_unif,
                         free_hi_res_unif,
                    continued_hi_res_unif,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)


io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'D'})

default_ds = 0.05 -- 0.02
num_continuate = 600
ds = default_ds

stepsize = 0.02
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
  io.write('        D : ', csurf.D, '\n')
  io.write('     ETA0 : ', csurf.ETA0, '\n')
  io.write('    U_MIN : ', csurf.U_MIN, '\n')
  io.write('    PHI_C : ', csurf.PHI_C, '\n')
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

  local theta = csurf.THETA_S
  local is_double = false
  for j = 2, #theta do
     if (theta[j] * theta[j-1] < 0) and (theta[j-1] > 0) then
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

    if nsurf.A > csurf.A then -- TEMPORARY HACK FOR RESTART
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
    freesurface, direction = htopography_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage('collect')

end

end -- if DO_US1

if DO_US2 then

--[[
Find the dip-like solutions? Then continue on to a pseudo-limiting configuration
with large aspect ratio of disturbance?
--]]

of = io.open(out_file_us2, "r")

if not of then
  io.write('Creating : ' .. out_file_us2 .. '.\n')
  io.flush()
  if not IN_SCRIPT_US then
    dofile(in_script_us)
    IN_SCRIPT_US = true
  end
  freesurface:create_file(out_file_us2)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_us2 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_us2)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_us2 .. '.\n')
    io.flush()
    if not IN_SCRIPT_US then
      dofile(in_script_us)
      IN_SCRIPT_US = true
    end
    freesurface:create_file(out_file_us2)
  else
    io.write('Loading : ' .. out_file_us2 .. '.\n')
    io.flush()
              hi_res_unif,
         free_hi_res_unif, 
    continued_hi_res_unif = htopography_load.load_last(out_file_us2,
                                                pert_to_unif_tables,
                                                           N_HI_RES)
    IN_SCRIPT_US = IN_SCRIPT_US or false
  end
end


freesurface = surface.new(    hi_res_unif,
                         free_hi_res_unif,
                    continued_hi_res_unif,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

default_ds = 0.05 -- 0.02
num_continuate = 400
ds = default_ds

stepsize = 0.02
max_stepsize = 2.5

freesurface:set_max_stepsize(max_stepsize)
freesurface:set_stepsize(stepsize)

direction = continuator.FORWARD

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

  if math.abs(csurf.A * csurf.L) > LARGEST_AREA then
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
    freesurface, direction = htopography_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage("collect")

end

end -- if DO_US2

if DO_DBL then

if not found_double then

of = io.open(out_file_dbl, "r")

if not of then
  io.write('Creating : ' .. out_file_dbl .. '.\n')
  io.flush()
  dofile(in_script_dbl)
  IN_SCRIPT_DBL = true
  freesurface:create_file(out_file_dbl)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_dbl ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_dbl)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_dbl .. '.\n')
    io.flush()
    dofile(in_script_dbl)
    IN_SCRIPT_DBL = true
    freesurface:create_file(out_file_dbl)
  else
    io.write('Loading : ' .. out_file_dbl .. '.\n')
    io.flush()
              hi_res_dbl,
         free_hi_res_dbl, 
    continued_hi_res_dbl = htopography_load.load_last(out_file_dbl,
                                             pert_to_double_tables,
                                                          N_HI_RES)
  end
end


freesurface = surface.new(     hi_res_dbl,
                          free_hi_res_dbl,
                     continued_hi_res_dbl,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)

--[[ Output file generation --]]

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

num_continuate = 400
default_ds = 0.05 -- 0.02
ds = default_ds

stepsize = 0.01 -- 0.02
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
  io.write('    PHI_C : ', csurf.PHI_C, '\n')
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

  if math.abs(csurf.A * csurf.L) > LARGEST_AREA then
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
    freesurface, direction = htopography_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage("collect")

end

end -- if not found_double

end -- if DO_DBL

if DO_DS2 then

of = io.open(out_file_ds2, "r")

if not of then
  io.write('Creating : ' .. out_file_ds2 .. '.\n')
  io.flush()
  if not IN_SCRIPT_DBL then
    dofile(in_script_dbl)
    IN_SCRIPT_DBL = true
  end
  freesurface:create_file(out_file_ds2)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_ds2 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_ds2)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_ds2 .. '.\n')
    io.flush()
    if not IN_SCRIPT_DBL then
      dofile(in_script_dbl)
      IN_SCRIPT_DBL = true
    end
    freesurface:create_file(out_file_ds2)
  else
    io.write('Loading : ' .. out_file_ds2 .. '.\n')
    io.flush()
              hi_res_dbl,
         free_hi_res_dbl, 
    continued_hi_res_dbl = htopography_load.load_last(out_file_ds2,
                                             pert_to_double_tables,
                                                          N_HI_RES)
    IN_SCRIPT_DS = false
  end
end

freesurface = surface.new(     hi_res_dbl,
                          free_hi_res_dbl,
                     continued_hi_res_dbl,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'FROUDE'})

num_continuate = 400
default_ds = 0.002 -- 0.015
ds = default_ds

stepsize = 0.002
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
  io.write('         A : ', csurf.A, '\n')
  io.write('    FROUDE : ', csurf.FROUDE, '\n')
  io.write('     U_MIN : ', csurf.U_MIN, '\n')
  io.write('     PHI_C : ', csurf.PHI_C, '\n')
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

  freesurface:append_file(out_file_ds2)

  if ((math.abs(csurf.A * csurf.L) > LARGEST_AREA) or
      (math.abs(csurf.A * csurf.L) < SMALLEST_AREA)) and
     (not initial_step) then
      mesgt.body = mesgt.body .. 
                   'Finished doube-hump solutions with reasonably large trench.\n' ..
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
    freesurface, direction = htopography_grid.regrid_hs(freesurface)
    local regrid_finish = os.clock()
    io.write('    Regrid time = ' .. tostring(regrid_finish - regrid_start) .. '\n')
    io.flush()
  end

  collectgarbage("collect")

end

end -- if DO_DS2


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
