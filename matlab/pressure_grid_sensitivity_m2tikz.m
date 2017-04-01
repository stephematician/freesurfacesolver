%% Consider the sensitivity of the solution to grid parameter CLUSTER_VAL
%
% Stephen Wade 01/07/2014

% So what this reveals is that, basically, the solutions converge very
% rapidly for small value of delta near 0.01 - 0.05. So we take delta
% = 0.02 and live with the results.

global PRINT_GRAPHS_TO_FILE;
global AXES_LIMITS;

AXES_LIMITS = [0, 0.09, 0, 2e-3]; % subject to change.

% Changed naming convention for files to improve sorting.
BATCHPLOT_STRETCH = 1.15;
% self explanatory
PRINT_GRAPHS_TO_FILE = true;
% This is the velocity at which 'stagnation' is assumed to occur
MIN_HORZ_VEL = 0.025;

disp('------------------------------ Process Gaussian Pressure : B = 2.8');
disp(' ');

prompt_str = input('Recompute results y/(n)? : ','s');
if strcmp(prompt_str, 'n') || strcmp(prompt_str, '')
  RECOMPUTE = false;
else
  RECOMPUTE = true;
end

if RECOMPUTE
  clear -regexp data_N
  RECOMPUTE = true;
end

if exist('profile_fig','var')
  if ishandle(profile_fig)
    close(profile_fig)
  end
  clear profile_fig
end

if exist('batch_fig_h','var')
  if ishandle(batch_fig_h)
    close(batch_fig_h)
  end
  clear batch_fig_h
end

if exist('branch_fig','var')
  if ishandle(branch_fig)
    close(branch_fig)
  end
  clear branch_fig
end

% Ask user for values of F to use
olddir = pwd();
try
  cd '../output/'
  disp('Available N numbers (allegedly):');
  a = dir('pressure_A010_GRID*.m');
  b = dir('pressure_A010_GRID*.mat');
  a1 = cell(1,length(a));
  for k=1:length(a)
    a1{k} = getfield(a, {k}, 'name');
  end
  b1 = cell(1,length(b));
  for k=1:length(b)
    b1{k} = getfield(b, {k}, 'name');
  end
  
  a2 = regexp(a1,'_GRID(?<number>\d+)\.', 'names');
  b2 = regexp(b1,'_GRID(?<number>\d+)\.', 'names');
  
  a3 = zeros(1,length(a2));
  for k=1:length(a2)
    a3(k) = str2double(a2{k}.number);
  end
  b3 = zeros(1,length(b2));
  for k=1:length(b2)
    b3(k) = str2double(b2{k}.number);
  end
  
  N_avail = unique([a3, b3]);
  N_avail = sort(N_avail)
  
catch err
  cd(olddir)
  throw(err)
end
cd(olddir)

Nindex = input('Input grid''s of interest BY INDEX, Nindex = ');

assert(isvector(Nindex), 'invalid selection for N')
assert(sum(ismember(Nindex,1:length(N_avail))) == length(Nindex), ...
       'invalid selection for N')

Nvec = N_avail(Nindex);

Nindex = input('Input ''exact'' grid of interest BY INDEX, Nindex = ');

assert(isscalar(Nindex), 'invalid selection for N')
assert(ismember(Nindex,1:length(N_avail)), ...
       'invalid selection for N')

Nexact = N_avail(Nindex);
set(0,'DefaultTextInterpreter','none');


%% Let's get to it.

disp('Ok, valid selection, let''s get to it')

batch_fig_h    = figure('units', 'centimeters', ...
                        'position', [0.5 9 7 5.25]);
                      
% use 'delta'.
f2 = @(in_d, i)((in_d.cluster_val{i} ./ in_d.phi_s_mid{i}(end)) *  ...
                filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL));

for N = Nvec
  old_wd = pwd();
  cd('../output/');
  try
    nstr = num2str(N);

    % Load data using some annoying code.
    
    evalstr = sprintf( ...
    ['if ~exist(''data_N%s'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''pressure_A010_GRID%s.m'', ''file'') || ...\n' ...
     '     exist(''pressure_A010_GRID%s.mat'', ''file'')\n' ...
     '    data_N%s = convert_and_load(''pressure_A010_GRID%s'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_N%s'', ''var'')\n' ...
     '  disp(''Check value of N (= %u) : '');\n' ...
     '  disp(length(data_N%s.x_s_mid{1})+1);\n' ...
     '  data_N%s.branch = ''UNF'';\n' ...
     'end\n'], ...
     nstr, nstr, nstr, nstr, nstr, nstr, N, nstr, nstr);
    eval(evalstr)
    
  evalstr = sprintf(...
  ['if exist(''data_N%s'', ''var'')\n' ...
   '  data_N%s.eta_s_pp = cell([1,length(data_N%s.A)]);\n' ...
   '  for k = 1:length(data_N%s.A)\n' ...
   '    data_N%s.eta_s_pp{k} = ...\n' ...
   '     spline(data_N%s.x_s_mid{k}, ...\n' ...
   '            data_N%s.eta_s_mid{k} + ...\n' ...
   '            (data_N%s.F{k}^2)*(1-(data_N%s.u_s_mid{k}(end)^2 + ...\n' ...
   '                                  data_N%s.v_s_mid{k}(end)^2))/2);\n' ...
   '  end\n' ...
   'end'], nstr, nstr, nstr, nstr, nstr, nstr, nstr, nstr, nstr, nstr);
    eval(evalstr)

  catch e
    cd(old_wd)
    throw(e)
  end
  
  cd(old_wd);
  
end

old_wd = pwd();
cd('../output/');
try
  nstr = num2str(Nexact);

  % Load data using some annoying code.

  evalstr = sprintf( ...
  ['if ~exist(''data_EXACT'', ''var'') || RECOMPUTE\n' ...
   '  if exist(''pressure_A010_GRID%s.m'', ''file'') || ...\n' ...
   '     exist(''pressure_A010_GRID%s.mat'', ''file'')\n' ...
   '    data_EXACT = convert_and_load(''pressure_A010_GRID%s'');\n' ...
   '  end\n' ...
   'end\n' ...
   'if exist(''data_EXACT'', ''var'')\n' ...
   '  disp(''Check value of N (= %u) : '');\n' ...
   '  disp(length(data_EXACT.x_s_mid{1})+1);\n' ...
   '  data_EXACT.branch = ''UNF'';\n' ...
   'end\n'], ...
   nstr, nstr, nstr, Nexact);
  eval(evalstr)

  evalstr = sprintf(...
  ['if exist(''data_EXACT'', ''var'')\n' ...
   '  data_EXACT.eta_s_pp = cell([1,length(data_EXACT.A)]);\n' ...
   '  for k = 1:length(data_EXACT.A)\n' ...
   '    data_EXACT.eta_s_pp{k} = ...\n' ...
   '     spline(data_EXACT.x_s_mid{k}, ...\n' ...
   '            data_EXACT.eta_s_mid{k} + ...\n' ...
   '            (data_EXACT.F{k}^2)*(1-(data_EXACT.u_s_mid{k}(end)^2 + ...\n' ...
   '                                    data_EXACT.v_s_mid{k}(end)^2))/2);\n' ...
   '  end\n' ...
   'end']);
  eval(evalstr)
catch e
  cd(old_wd)
  throw(e)
end

cd(old_wd);

j = 1;

c = blue_black(length(Nvec));


for N = Nvec
  this_color = c(j,:);

  figure(batch_fig_h)
  hold on
  nstr = num2str(N);
  
  evalstr = sprintf( ...
  ['if exist(''data_N%s'',''var'')\n' ...
   '  [~] = batch_plot_m2tikz(data_N%s, ...\n' ...
   '             f2, @(x,y)(unforced_integral_check_b(x,y,data_EXACT)), ...\n' ...
   '             ''-'', ''color'', this_color);\n' ...
   'end'], ...
  nstr, nstr);

  eval(evalstr)
  hold off
  
  j = j + 1;
  
end

axis(AXES_LIMITS)


ylabel('$E$', 'rotation', 0);
xlabel('$\delta_1$');

grid on
box on

if PRINT_GRAPHS_TO_FILE
  old_wd = pwd;
  cd('../output/')
  file_name = 'grid_sensitivity_summary.tikz';
  matlab2tikz(file_name, 'width','\figurewidth', ...
              'parseStrings', false, ...
              'figurehandle', batch_fig_h);
  cd(old_wd); 
end