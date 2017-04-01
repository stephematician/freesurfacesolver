--[[
Compute continuation curve of pressure disturbance with fixed B, and fixed
Froude number.

Performs two or three continuations. The first is from uniform stream, and it
continues clockwise around until the velocity at the crest gets too small, or
it hits double hump solutions with a small A. It then performs a
continuation from the original uniform stream to the 'dip' like solutions.
If no double hump solutions were found in the first continuation, then this
finds a double hump solution and then performs a clockwise continuation.

rev03 - homotopy available, symmetric formulation, with dbeta dependent
        upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 22/05/2012
Modified - 22/05/2013

History:
22/5/2013 Created file, adapted from gtopography_F110_B500.lua

]]--

require "surface"

EMAILME = true
OVERWRITE = true
DO_STOKES = true
DO_ALT = false

in_script  = 'input/pressure_A000_B280_init.lua'

out_file_stokes = 'output/pressure_A000_STO_B280.m'
out_file_alt = 'output/pressure_A000_ALT_B280.m'

--[[ This file creates the initial_surf, is_free_surf, etc tables. ]]--
dofile(in_script)

--[[ Email notification sent when computation is completed. ]]--
if EMAILME then
  require "luarocks.loader"
  smtp = require("socket.smtp")

  from = "<stephen.wade@adelaide.edu.au>"
  smtpserver = "smtp.adelaide.edu.au"

  rcpt = {
    "<stephematician@gmail.com>"
  }

  mesgt = {
    headers = {
      to = "Stephen Wade <stephematician@gmail.com>",
      from = "Stephen Wade <stephen.wade@adelaide.edu.au>",
      subject = "pressure calculation - stokes" 
    },
    body = "This process should be complete.\n\n" ..
           "Notes: homotopy available, symmetric formulation, with beta and " ..
           "dbeta dependent upon value of u.\n"
  }
end

if DO_STOKES then

freesurface = surface.new(   initial_surf,
                             is_free_surf,
                        is_continued_surf,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)


--[[ Output file generation ]]--

of = io.open(out_file_stokes, "r")

if not of then
   freesurface:create_file(out_file_stokes)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_stokes)
    freesurface:create_file(out_file_stokes) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'U_MIN', 'ETA0'})

num_continuate = 400
ds = 0.008

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
  io.write('            A : ', csurf.A, '\n')
  io.write('       froude : ', csurf.FROUDE, '\n')
  io.write('         eta0 : ', csurf.ETA0, '\n')
  io.write('         umin : ', csurf.U_MIN, '\n')
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
    if EMAILME then
      mesgt.body = mesgt.body .. 
                   'Finished due to error in continuation method.\n'
    end
    error('Error detected in continuation.')
    break
  end

  freesurface:append_file(out_file_stokes)

  if csurf.U_MIN < 0.01 then
    if EMAILME then
      mesgt.body = mesgt.body .. 
                   'Finished first continuation branch at a Stokes configuration.\n' ..
                   'Took ' .. i .. ' iterations.\n'
    end

    break
  end

  if csurf.U_MIN < 0.15 then
    ds = 0.0025
  elseif csurf.U_MIN < 0.40 then
    ds = 0.0050
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

    if nsurf.U_MIN > csurf.U_MIN then
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

end -- DO_STOKES

--[[

]]--

if DO_ALT then

freesurface = surface.new(   initial_surf,
                             is_free_surf,
                        is_continued_surf,
           pressure_residuals.residual_hs,
     pressure_residuals.compute_output_hs)

of = io.open(out_file_alt, "r")

if not of then
   freesurface:create_file(out_file_alt)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file_alt)
    freesurface:create_file(out_file_alt) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'U_MIN', 'ETA0'})

num_continuate = 400
ds = 0.025

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
  io.write('            A : ', csurf.A, '\n')
  io.write('       froude : ', csurf.FROUDE, '\n')
  io.write('         eta0 : ', csurf.ETA0, '\n')
  io.write('         umin : ', csurf.U_MIN, '\n')
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
    if EMAILME then
      mesgt.body = mesgt.body .. 
                   'Finished due to error in continuation method.\n'
    end
    error('Error detected in continuation.')
    break
  end

  freesurface:append_file(out_file_alt)

  if math.abs(csurf.A * csurf.B) > 15 then
    if EMAILME then
      mesgt.body = mesgt.body .. 
                   'Finished dip-like solutions with reasonably large dip.\n' ..
                   'Took ' .. i .. ' iterations.\n'
    end

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
    freesurface, direction = pressure_grid.regrid_hs(freesurface)
  end

end

end -- if DO_ALT

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
