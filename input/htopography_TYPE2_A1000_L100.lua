--[[
Compute continuation curve of topography disturbance with fixed L, B and A.

Box is defined via the use of an approximation to a step function given by
the hyperbolic tangent.

Performs a single continuation to find a limiting solution with a single crest
with A fixed.

Author   - Stephen Wade
Created  - 26/08/2014
Modified - 17/01/2015

History:
22/05/2013  Created file, adapted from gtopography_F110_B500.lua
14/12/2013  Adapted from pressure_F127_B280.lua
26/08/2014  Adapted from pressure_TYPE4_A1000_B280.lua
07/09/2014  Adapted from htopography_TYPE4_A1000_B280.lua
16/01/2015  Adapted from htopography_U004_L100.lua
]]--

require "surface"
require "htopography_residuals"
require "htopography_grid"
require "htopography_load"

dofile('input/htopography_find_type2_solitary.lua')

EMAILME = false
OVERWRITE = false

N_HI_RES = 1000 -- 1400
N_LO_RES = 240 -- 600 should do it.

CLUSTER_DPHI_VAL = 0.01

-- Domain and topography
L       = 1.00
B       = 10.00 -- ?
PHI_LIM = 10

-- Starting and finishing
TYPE2_A = 0.1000 -- ?
TYPE2_FROUDE = 1.22
STOKES_U_MAX = 0.5
STOKES_U_MIN = 0.04

round = function(x)
  if x%2 ~= 0.5 then
    return math.floor(x+0.5)
  end
  return x-0.5
end

ASTR = string.format('%03u',round(math.abs(TYPE2_A*100)))
USTR = string.format('%03u',round(STOKES_U_MIN*100))

in_script  = 'input/htopography_type2_solitary_init.lua' -- alt?

out_file = 'output/htopography_A' .. ASTR .. '_U' .. USTR .. '_TYPE2.m'

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
    subject = 'htopography calculation htopography_A' .. ASTR .. '_U' .. USTR .. '_TYPE4.m' 
  },
  body = "This process should be complete.\n\n" ..
         "Notes: The usual stuff, all as per thesis.\n"
}

--[[ Output file generation ]]--

of = io.open(out_file, "r")

if not of then
  io.write('Creating : ' .. out_file .. '.\n')
  io.flush()
  dofile(in_script)
  freesurface:create_file(out_file)
else
  if OVERWRITE then
    io.write('Deleting : ' .. out_file ..' ...')
    io.flush()
    of:close()
    os.remove(out_file)
    io.write('done\n')
    io.flush()
    io.write('Creating : ' .. out_file .. '.\n')
    io.flush()
    dofile(in_script)
    freesurface:create_file(out_file)
  else
    io.write('Loading : ' .. out_file .. '.\n')
    io.flush()
              hi_res_sol,
         free_hi_res_sol, 
    continued_hi_res_sol = htopography_load.load_last(out_file,
                                           type2_single_tables,
                                                      N_HI_RES)
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

freesurface:set_length_vars({'U_MIN', 'ETA0'})

num_continuate = 400
ds = 0.015

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
  io.write('     U_MIN : ', csurf.U_MIN, '\n')
  io.write('      ETA0 : ', csurf.ETA0, '\n')
  io.write('    FROUDE : ', csurf.FROUDE, '\n')
  io.write('         A : ', csurf.A, '\n')

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
    freesurface, direction = htopography_grid.regrid_hs(freesurface)
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
