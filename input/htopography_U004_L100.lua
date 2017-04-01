--[[
Compute continuation curve of smooth box shaped disturbance with fixed width,
and fixed U at the crest.

Box is defined via the use of an approximation to a step function given by
the hyperbolic tangent.

Performs three or four continuations. The first is from a solitary wave, with 
small u at the crest and a given Froude number. It traverses to increasing
Froude number and ends once it hits a large enough topographical disturbance.
The second is from the same limiting solitary wave, but it traverses along
decreasing Froude number until again a large topographical disturbance is found.
The last two are similar, but start from a double hump solitary wave with a
small u at the crest. 


Author   - Stephen Wade
Created  - 07/12/2012
Modified - 26/07/2014

History:
26/07/2014 Adapted from htopography_F128_L200.lua
...
07/12/2012 Created original script.
--]]

require "surface"
require "htopography_residuals"
require "htopography_grid"
require "htopography_load"

dofile('input/htopography_find_stokes_solitary.lua')
dofile('input/htopography_find_stokes_double.lua')

DO_SS1 = false
DO_SS2 = false
DO_DS1 = true
DO_DS2 = false

EMAILME = false
OVERWRITE = false

-- Grid
N_HI_RES = 1800
N_LO_RES = 480

CLUSTER_DPHI_VAL = 0.01

-- Domain and topography
L            = 1.00
B            = 10.00 -- ?
FROUDE       = 1.28
DOUBLE_A     = -0.02
PHI_LIM      = 10
STOKES_U_MIN = 0.04

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

USTR = string.format('%03u',round(STOKES_U_MIN*100))
LSTR = string.format('%03u',round(L*100))

-- Starting and finishing
UNIFORM_A     = -0.04
LIMIT_DBL_A   = 0.0005*L
LARGEST_AREA  = 5.0
SMALLEST_AREA = 0.0025
LARGEST_ETA0  = 2.5

in_script_ss = 'input/htopography_stokes_solitary_init.lua'
in_script_ds = 'input/htopography_stokes_double_init.lua'


out_file_ss1 = 'output/htopography_U' .. USTR .. '_SS1_L' .. LSTR .. '.m'
out_file_ss2 = 'output/htopography_U' .. USTR .. '_SS2_L' .. LSTR .. '.m'
out_file_ds1 = 'output/htopography_U' .. USTR .. '_DS1_L' .. LSTR .. '.m'
out_file_ds2 = 'output/htopography_U' .. USTR .. '_DS2_L' .. LSTR .. '.m'

--[[ Email notification sent when computation is completed. --]]
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
    subject = "htopography calculation" 
  },
  body = "Body:\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}


if DO_SS1 then

of = io.open(out_file_ss1, "r")

if not of then
  io.write('Creating : ' .. out_file_ss1 .. '.\n')
  io.flush()
  dofile(in_script_ss)
  IN_SCRIPT_SS = true
  freesurface:create_file(out_file_ss1)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_ss1 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_ss1)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_ss1 .. '.\n')
    io.flush()
    dofile(in_script_ss)
    freesurface:create_file(out_file_ss1)
    IN_SCRIPT_SS = true
  else
    io.write('Loading : ' .. out_file_ss1 .. '.\n')
    io.flush()
              hi_res_sol,
         free_hi_res_sol, 
    continued_hi_res_sol = htopography_load.load_last(out_file_ss1,
                                              stokes_single_tables,
                                                          N_HI_RES)
    IN_SCRIPT_SS = false
  end
end


freesurface = surface.new(     hi_res_sol,
                          free_hi_res_sol,
                     continued_hi_res_sol,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)


io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'FROUDE'})

default_ds = 0.04 -- 0.05
num_continuate = 600
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
  io.write('         A : ', csurf.A, '\n')
  io.write('      ETA0 : ', csurf.ETA0, '\n')
  io.write('     U_MIN : ', csurf.U_MIN, '\n')
  io.write('    FROUDE : ', csurf.FROUDE, '\n')
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

  freesurface:append_file(out_file_ss1)

  if (math.abs(csurf.A * csurf.L) > LARGEST_AREA) or
     (csurf.ETA0 > LARGEST_ETA0) then
    mesgt.body = mesgt.body .. 
                 'Finished first continuation branch at a Stokes configuration.\n' ..
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

  collectgarbage('collect')

end

end -- if DO_SS1

if DO_SS2 then

of = io.open(out_file_ss2, "r")

if not of then
  io.write('Creating : ' .. out_file_ss2 .. '.\n')
  io.flush()
  if not IN_SCRIPT_SS then
    dofile(in_script_ss)
    IN_SCRIPT_SS = true
  end
  freesurface:create_file(out_file_ss2)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_ss2 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_ss2)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_ss2 .. '.\n')
    io.flush()
    if not IN_SCRIPT_SS then
      dofile(in_script_ss)
      IN_SCRIPT_SS = true
    end
    freesurface:create_file(out_file_ss2)
  else
    io.write('Loading : ' .. out_file_ss2 .. '.\n')
    io.flush()
              hi_res_sol,
         free_hi_res_sol, 
    continued_hi_res_sol = htopography_load.load_last(out_file_ss2,
                                              stokes_single_tables,
                                                          N_HI_RES)
    IN_SCRIPT_SS = IN_SCRIPT_SS or false
  end
end


freesurface = surface.new(     hi_res_sol,
                          free_hi_res_sol,
                     continued_hi_res_sol,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)


--[[ Output file generation --]]

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'FROUDE'})

default_ds = 0.04 -- 0.05
num_continuate = 600
ds = default_ds

stepsize = 0.01
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
  io.write('         A : ', csurf.A, '\n')
  io.write('      ETA0 : ', csurf.ETA0, '\n')
  io.write('     U_MIN : ', csurf.U_MIN, '\n')
  io.write('    FROUDE : ', csurf.FROUDE, '\n')
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

  freesurface:append_file(out_file_ss2)

  if (math.abs(csurf.A * csurf.L) > LARGEST_AREA) or
     (csurf.ETA0 > LARGEST_ETA0) then
    mesgt.body = mesgt.body .. 
                 'Finished first continuation branch at a Stokes configuration.\n' ..
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

  collectgarbage('collect')

end

end -- if DO_SS2

--print(mesgt.headers.subject)
--print(mesgt.body)

io.flush()

if DO_DS1 then

of = io.open(out_file_ds1, "r")

if not of then
  io.write('Creating : ' .. out_file_ds1 .. '.\n')
  io.flush()
  dofile(in_script_ds)
  IN_SCRIPT_DS = true
  freesurface:create_file(out_file_ds1)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file_ds1 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_ds1)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_ds1 .. '.\n')
    io.flush()
    dofile(in_script_ds)
    freesurface:create_file(out_file_ds1)
    IN_SCRIPT_DS = true
  else
    io.write('Loading : ' .. out_file_ds1 .. '.\n')
    io.flush()
              hi_res_dbl,
         free_hi_res_dbl, 
    continued_hi_res_dbl = htopography_load.load_last(out_file_ds1,
                                              stokes_double_tables,
                                                          N_HI_RES)
  -- TEMPORARY CHANGES
freesurface = surface.new(     hi_res_dbl,
                          free_hi_res_dbl,
                     continued_hi_res_dbl,
        htopography_residuals.residual_hs,
  htopography_residuals.compute_output_hs)
    io.write('Deleting : ' .. out_file_ds1 ..' ...')
    io.flush()
    of:close()
    os.remove(out_file_ds1)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file_ds1 .. '.\n')
    io.flush()
    freesurface:create_file(out_file_ds1)
    -- TEMPORARY CHANGES
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
default_ds = 0.015
ds = default_ds

stepsize = 2
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

  freesurface:append_file(out_file_ds1)

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

end -- if DO_DS1

if DO_DS2 then

of = io.open(out_file_ds2, "r")

if not of then
  io.write('Creating : ' .. out_file_ds2 .. '.\n')
  io.flush()
  if not IN_SCRIPT_DS then
    dofile(in_script_ds)
    IN_SCRIPT_DS = true
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
    if not IN_SCRIPT_DS then
      dofile(in_script_ds)
      IN_SCRIPT_DS = true
    end
    freesurface:create_file(out_file_ds2)
  else
    io.write('Loading : ' .. out_file_ds2 .. '.\n')
    io.flush()
              hi_res_dbl,
         free_hi_res_dbl, 
    continued_hi_res_dbl = htopography_load.load_last(out_file_ds2,
                                              stokes_double_tables,
                                                          N_HI_RES)
    IN_SCRIPT_DS = IN_SCRIPT_DS or false
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

freesurface:set_length_vars({'A', 'FROUDE'})

num_continuate = 400
default_ds = 0.025
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
