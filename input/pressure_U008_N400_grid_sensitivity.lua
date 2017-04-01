--[[
Testing the sensitivity of the solution to the grid clustering values
that I can choose.

Computes for varying CLUSTER_VALS the nearly unforced almost-highest
solitary wave

Performs one continuation from a perturbation to a solitary wave with U=0.008
with F = 1.29089, from CLUSTER_VAL = 0.005 to ???

rev03 - homotopy available, symmetric formulation, with dbeta dependent
        upon value of u^2, and beta dependent on 1/u (averaged).

Author   - Stephen Wade
Created  - 23/06/2014
Modified - 

History:
20/06/2014 Adapted from pressure_U010_B280.lua

--]]

require "surface"
require "pressure_residuals"
require "pressure_grid"

EMAILME = false
OVERWRITE = true

N_HI_RES = 400 -- 800
N_LO_RES = 200 -- 400

B                = 2.80
FROUDE           = 1.29089
PHI_LIM          = 11
CLUSTER_DPHI_VAL = 0.005
A                = 0.00001

STOKES_U_MIN = 0.08
MAX_CLUSTER_VAL = 2

display_initial = false

in_script  = 'input/pressure_stokes_solitary_init.lua'

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

ASTR = string.format('%3u',round(A*1E6))

out_file = 'output/pressure_A' .. ASTR ..  '_GRID' .. tostring(N_HI_RES) .. '.m'

-- This file creates the initial_surf, is_free_surf, etc tables.

dofile(in_script)

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
  body = "This process should be complete.\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}

freesurface = surface.new(    hi_res_sol,
                         free_hi_res_sol,
                    continued_hi_res_sol,
          pressure_residuals.residual_hs,
    pressure_residuals.compute_output_hs)

-- Output file generation

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

freesurface:set_length_vars({'U_MIN', 'CLUSTER_VAL'})

num_continuate = 400
ds = 0.03

stepsize = 0.0001
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
  io.write('    CLUSTER_VAL : ', csurf.CLUSTER_VAL, '\n')
  io.write('           ETA0 : ', csurf.ETA0, '\n')
  io.write('          U_MIN : ', csurf.U_MIN, '\n')
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

  if (math.abs(csurf.CLUSTER_VAL) > MAX_CLUSTER_VAL) then
    mesgt.body = mesgt.body .. 
                 '\nFinished continuation branch.\n' ..
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

    if nsurf.CLUSTER_VAL < csurf.CLUSTER_VAL then
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

io.write('Finished continuation.\n')
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
