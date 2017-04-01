--[[
Compute continuation curve of pressure disturbance with fixed B, and fixed
u at crest.

Performs three continuations. The first is from a perturbation to a solitary
wave with F=1.2909 with small U at the crest, and then it finds the positively
forced solutions with fixed u up until it hits a large dip/bump. Then it finds
all the negatively forced solutions from the same initial profile as the first
process. Then, finally, it finds the double hump branch, starting with a double
hump solution with small u at the crest and small -ve A, and then follows
this branch around to a reasonably large dip/bump.

rev03 - homotopy available, symmetric formulation, with dbeta dependent
        upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 22/05/2013
Modified - 14/12/2013

History:
22/5/2013  Created file, adapted from gtopography_F110_B500.lua
14/12/2013 Adapted from pressure_F127_B280.lua

]]--

require "surface"
require "pressure_residuals"
require "pressure_grid"

EMAILME = false
OVERWRITE = true

DO_SL1 = false
DO_SL2 = false
DO_DBL = true

N = 800       -- 800
N_ROUGH = 400 -- 400

B                = 2.80
FROUDE           = 1.29089
PHI_LIM          = 11
CLUSTER_DPHI_VAL = 0.0025

HELP_BIFUR_A     = -0.25
MAX_AREA         = 0.5*B

STOKES_U_MIN = 0.08

display_initial = false

in_script  = 'input/pressure_stokes_solitary_init.lua'
in_script_dbl  = 'input/pressure_stokes_double_init.lua'

out_file_sl1 = 'output/pressure_U008_SL1_B280.m'
out_file_sl2 = 'output/pressure_U008_SL2_B280.m'
out_file_dbl = 'output/pressure_U008_DBL_B280.m'

--[[ This file creates the initial_surf, is_free_surf, etc tables. ]]--
if DO_SL1 or DO_SL2 then
  dofile(in_script)
end
--[[ Email notification sent when computation is completed. ]]--
if EMAILME then
  require "luarocks.loader"
  smtp = require("socket.smtp")
end

from = "<stephen.wade@adelaide.edu.au>"
smtpserver = "smtp.adelaide.edu.au"

rcpt = {
  "<stephematician@gmail.com>"
}

mesgt = {
  headers = {
    to = "Stephen Wade <stephematician@gmail.com>",
    from = "Stephen Wade <stephen.wade@adelaide.edu.au>",
    subject = "gtopography calculation" 
  },
  body = "This process should be complete.\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}

if DO_SL1 then

freesurface = surface.new(    initial_sol,
                              is_free_sol,
                         is_continued_sol,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)


--[[ Output file generation ]]--

of = io.open(out_file_sl1, "r")

if not of then
   freesurface:create_file(out_file_sl1)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_sl1)
    freesurface:create_file(out_file_sl1) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

num_continuate = 400
ds = 0.01

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
  io.write('          a : ', csurf.A, '\n')
  io.write('       eta0 : ', csurf.ETA0, '\n')
  io.write('          F : ', csurf.FROUDE, '\n')
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

  freesurface:append_file(out_file_sl1)

  if (math.abs(csurf.A * csurf.B) > MAX_AREA) then
    mesgt.body = mesgt.body .. 
                 '\nFinished first continuation branch at a large bump.\n' ..
                 'Took ' .. i .. ' iterations.\n'
    break
  end


  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  io.write(msg)
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
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
  end

end

end -- DO_SL1

--[[
Find the dip-like solutions? Then continue on to a pseudo-limiting configuration
with large aspect ratio of disturbance?
]]--

if DO_SL2 then

freesurface = surface.new(    initial_sol,
                              is_free_sol,
                         is_continued_sol,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

of = io.open(out_file_sl2, "r")

if not of then
   freesurface:create_file(out_file_sl2)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_sl2)
    freesurface:create_file(out_file_sl2) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'A', 'ETA0'})

num_continuate = 400
ds = 0.005

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
found_double = false

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
  io.write('          a : ', csurf.A, '\n')
  io.write('       eta0 : ', csurf.ETA0, '\n')
  io.write('          F : ', csurf.FROUDE, '\n')
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

  freesurface:append_file(out_file_sl2)

  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  io.write(msg)
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
    local nsurf = surface.unpack_vec(freesurface.unknowns,
                                      freesurface.initial,
                                      freesurface.is_free,
                                 freesurface.is_continued,
                                    freesurface.ind_table)

    if nsurf.FROUDE < csurf.FROUDE then
      mesgt.body = mesgt.body .. 
                   '\nFinished dip-like solutions with reasonably large dip.\n' ..
                   'Took ' .. i .. ' iterations.\n'
      break
    end
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
  end

end

end -- if DO_SL2

if DO_DBL then

dofile(in_script_dbl)

freesurface = surface.new(    initial_dbl,
                              is_free_dbl,
                         is_continued_dbl,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

--[[ Output file generation ]]--

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
ds = 0.005

stepsize = 0.01
max_stepsize = 4.0

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
  io.write('          A : ', csurf.A, '\n')
  io.write('       ETA0 : ', csurf.ETA0, '\n')
  io.write('      U_MIN : ', csurf.U_MIN, '\n')
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

  res, msg = freesurface:progress_interp(ds, direction, freesurface)
  io.write(msg)
  io.flush()

  local nsurf = surface.unpack_vec(freesurface.unknowns,
                                    freesurface.initial,
                                    freesurface.is_free,
                               freesurface.is_continued,
                                  freesurface.ind_table)
  if (not initial_step) then
     if ((nsurf.A > csurf.A) and csurf.A < HELP_BIFUR_A) or
      (math.abs(csurf.A * csurf.B) > MAX_AREA) then
      mesgt.body = mesgt.body .. 
                   'Finished doube-hump solutions with either reasonably large trench ' ..
                   'or the solutions started to double back on themselves.\n' ..
	           'Took ' .. (i+1) .. ' iterations.\n'
      break
   end

   freesurface, direction = pressure_grid.regrid_hs(freesurface)

  else
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
  end

end

end -- DO DBL

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
