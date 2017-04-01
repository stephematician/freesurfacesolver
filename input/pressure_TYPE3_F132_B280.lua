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

N_HI_RES = 2400 -- 800
N_LO_RES = 600 -- 380 should do it.

B                = 2.80
FROUDE           = 1.32
PHI_LIM          = 11
CLUSTER_DPHI_VAL = 0.005

display_initial = false

-- Starting and finishing
STOKES_U_MAX = 0.33
STOKES_U_MIN = 0.015

HELPER_A = -0.15668

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

FSTR = string.format('%3u',round(FROUDE*100))
USTR = string.format('%03u',round(STOKES_U_MIN*100))

in_script  = 'input/pressure_stokes_type3_init.lua'

out_file = 'output/pressure_F' .. FSTR .. '_U' .. USTR .. '_TYPE3.m'

--[[ This file creates the initial_surf, is_free_surf, etc tables. ]]--

dofile(in_script)

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
    subject = 'pressure calculation pressure_F' .. FSTR .. '_U' .. USTR .. '_TYPE3.m' 
  },
  body = "This process should be complete.\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}

free_hi_res_sol.FROUDE = true
free_hi_res_sol.A = false

continued_hi_res_sol.A = false
continued_hi_res_sol.U_MIN = true

freesurface = surface.new(    hi_res_sol,
                         free_hi_res_sol,
                    continued_hi_res_sol,
          pressure_residuals.residual_hs,
    pressure_residuals.compute_output_hs)

--[[ Output file generation ]]--

of = io.open(out_file, "r")

if not of then
   freesurface:create_file(out_file)
else
  if OVERWRITE then
    of:close()
    os.remove(out_file)
    freesurface:create_file(out_file) 
  end
end

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

freesurface:set_length_vars({'U_MIN', 'ETA0'})

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
  io.write('     U_MIN : ', csurf.U_MIN, '\n')
  io.write('      ETA0 : ', csurf.ETA0, '\n')
  io.write('    FROUDE : ', csurf.FROUDE, '\n')
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

  freesurface:append_file(out_file)

  if (math.abs(csurf.U_MIN) < STOKES_U_MIN) then
    mesgt.body = mesgt.body .. 
                 '\nFinished first continuation branch with non-stokes crest.\n' ..
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
