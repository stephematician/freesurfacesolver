--[[
Testing the sensitivity of the solution to the domain size
that I can choose.

Computes for varying PHI_M the nearly unforced almost-highest
solitary wave

Performs one continuation from a perturbation to a solitary wave with U=0.005
with F = 1.29089, from PHI_M = 7 to 12

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

N_REF_RES = 548 -- 800
N_LO_RES = 200 -- 400

B                = 2.80
FROUDE           = 1.32
PHI_REF          = 11
PHI_MIN          = 7
PHI_MAX          = 13
CLUSTER_DPHI_VAL = 0.01
HELPER_A         = -0.15668

TOTAL_RUN = 15

-- Starting and finishing
STOKES_U_MAX = 0.05
STOKES_U_MIN = 0.05
MAX_CLUSTER_VAL = 2

display_initial = false

in_script  = 'input/pressure_stokes_type3_init.lua'

out_file = 'output/pressure_TYPE3_PHIM_N' .. tostring(N_REF_RES) .. '.m'

N_MIN = math.floor(PHI_MIN * N_REF_RES / PHI_REF)
N_MAX = math.ceil(PHI_MAX * N_REF_RES / PHI_REF)
N_STEP = math.ceil((N_MAX - N_MIN) / TOTAL_RUN)

N_HI_RES = N_MIN
PHI_LIM = N_MIN * PHI_REF / N_REF_RES

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
    subject = "pressure phi_m sensitivity calculation" 
  },
  body = "This process should be complete.\n\n" ..
         "Notes: homotopy available, symmetric formulation, with beta and " ..
         "dbeta dependent upon value of u.\n"
}

continued_hi_res_sol.A = false
continued_hi_res_sol.U_MIN = true

hi_res_sol.U_MIN  = STOKES_U_MIN
hi_res_sol.FROUDE = FROUDE

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

freesurface:newton(freesurface)

freesurface:append_file(out_file)

io.write('---------------------\n')
io.write('Starting continuation\n')
io.flush()

-- Want to loop from N = N_MIN to N_MAX

for N = (N_MIN+N_STEP), (N_MIN+N_STEP*TOTAL_RUN), N_STEP do
  local outstr = 'Current N : ' .. tostring(N-1)
  local filstr = '-'
  io.write(outstr .. '\n')
  io.write(filstr:rep(outstr:len()) .. '\n')
    

  local new_res_sol = surface.unpack_vec(freesurface.unknowns,
                                          freesurface.initial,
                                          freesurface.is_free,
                                     freesurface.is_continued,
                                        freesurface.ind_table)
  local free_new_res_sol      = pressure_grid.deepcopy(freesurface.is_free)
  local continued_new_res_sol = pressure_grid.deepcopy(freesurface.is_continued)
  
  io.write('  Parameters for previous solution.\n')
  io.write('       N :')
  for k,v in pairs(new_res_sol.N_SUB) do
     io.write(' ' .. tostring(v))
  end
  io.write(',\n')
  io.write('       A : ', new_res_sol.A, ',\n')
  io.write('       B : ', new_res_sol.B, ',\n')
  io.write('       F : ', new_res_sol.FROUDE, ',\n')
  io.write('  LAMBDA : ', new_res_sol.LAMBDA, ',\n')
  io.write('    ETA0 : ', new_res_sol.ETA0, ',\n')
  io.write('   U_MIN : ', new_res_sol.U_MIN, '.\n\n')
  
  new_res_sol.BETA_SUB[#new_res_sol.BETA_SUB] = N
  new_res_sol.PHI_SUB[#new_res_sol.PHI_SUB] = N * PHI_REF / N_REF_RES
  new_res_sol.N_SUB[#new_res_sol.N_SUB] = new_res_sol.N_SUB[#new_res_sol.N_SUB] + N_STEP
  -- hopefully dphi sub comes from the solution.
                                               
  local y01 = new_res_sol.V_S[N-N_STEP]
  local y02 = new_res_sol.V_S[N-(N_STEP+1)]
  local l0 = math.log(y01/y02)
  local A0 = y01 * math.exp(-l0 * (N-N_STEP))
  
  for k = 1, N_STEP do
    new_res_sol.V_S[N-N_STEP+k] = A0 * math.exp(l0 * (N-N_STEP+k))
  end

  free_new_res_sol.V_S = {}
  for k = 1, N do
    free_new_res_sol.V_S[k] = true
  end
  free_new_res_sol.V_S[1] = false

  continued_new_res_sol.THETA_S = {}
  for k = 1, N do
    continued_new_res_sol.THETA_S[k] = false
  end
  
  freesurface = surface.new(    new_res_sol,
                           free_new_res_sol,
                      continued_new_res_sol,
             pressure_residuals.residual_hs,
       pressure_residuals.compute_output_hs)
       
  res = freesurface:newton(freesurface)
  
  freesurface:append_file(out_file)

  if res == continuator.ERROR then
    error('Unexpected error detected in continuator:progress().')
    break
  end

  if res == continuator.DIVERGED then
    error('Unexpected divergence detected in continuator:progress().')
    break
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
