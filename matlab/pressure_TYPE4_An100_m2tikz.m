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

disp('---------------------- Process Gaussian Pressure : TYPE 4 Solution');
disp(' ');

prompt_str = input('Recompute results y/(n)? : ','s');
if strcmp(prompt_str, 'n') || strcmp(prompt_str, '')
  RECOMPUTE = false;
else
  RECOMPUTE = true;
end

if RECOMPUTE
  clear -regexp cb_xd_|cb_yd_|data_F|xd_|yd_|PLOT_F
  RECOMPUTE = true;
end

if exist('profile_fig','var')
  if ishandle(profile_fig)
    close(profile_fig)
  end
  clear profile_fig
end

if exist('theta_fig','var')
  if ishandle(theta_fig)
    close(theta_fig)
  end
  clear theta_fig
end

if exist('batch_fig_h','var')
  if ishandle(batch_fig_h)
    close(batch_fig_h)
  end
  clear batch_fig_h
end

if exist('branch_fig_E','var')
  if ishandle(branch_fig_E)
    close(branch_fig_E)
  end
  clear branch_fig_E
end

if exist('branch_fig_h','var')
  if ishandle(branch_fig_h)
    close(branch_fig_h)
  end
  clear branch_fig_h
end

% self explanatory
PRINT_GRAPHS_TO_FILE = true;

% This is the velocity at which 'stagnation' is assumed to occur
MIN_HORZ_VEL = 0.001; %0.1899

date_str = strrep(date(),'-','_');

% Ask user for values of F to use
olddir = pwd();
try
  cd '../output/'
  disp('Available N numbers (allegedly):');
  a = dir('pressure_TYPE4_A*.m');
  b = dir('pressure_TYPE4_A*.mat');
  a1 = cell(1,length(a));
  for k=1:length(a)
    a1{k} = getfield(a, {k}, 'name');
  end
  b1 = cell(1,length(b));
  for k=1:length(b)
    b1{k} = getfield(b, {k}, 'name');
  end
  
  a2 = regexp(a1,'_A(?<number>\d+)\.', 'names');
  b2 = regexp(b1,'_A(?<number>\d+)\.', 'names');
  
  a3 = zeros(1,length(a2));
  for k=1:length(a2)
    a3(k) = str2double(a2{k}.number);
  end
  b3 = zeros(1,length(b2));
  for k=1:length(b2)
    b3(k) = str2double(b2{k}.number);
  end
  
  A_avail = unique([a3, b3]);
  A_avail = sort(N_avail)
  
catch err
  cd(olddir)
  throw(err)
end
cd(olddir)

Aindex = input('Input A of interest BY INDEX : ');

assert(isvector(Aindex), 'invalid selection for A')
assert(sum(ismember(Aindex,1:length(A_avail))) == length(Aindex), ...
       'invalid selection for A')

Avec = A_avail(Aindex);

set(0,'DefaultTextInterpreter','none');

%% Let's get to it.

disp('Ok, valid selection, let''s get to it')

batch_fig_h    = figure('units', 'centimeters', ...
                        'position', [0.5 9 7 5.25]);
                      
profile_fig = figure('units', 'centimeters', ...
                     'position', [10.5 9 7 5.25]);

branch_fig = figure('units', 'centimeters', ...
                    'position', [20.5 9 7 5.25]);

% f1 = omega
f1 = @(in_d, i)((1 - (in_d.F{i} * min(in_d.u_s_mid{i}))^2) * ...
                filter_trick((min(in_d.u_s_mid{i})) > MIN_HORZ_VEL));
% f2 = F
f2 = @(in_d, i)((in_d.F{i}) *  ...
                filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL));
% f3 = corrected eta0
f3 = @(in_d, i)((in_d.eta_s_mid{i}(1) + ...
                 (data_EXACT.F{k}^2) * ...
                 (1-(data_EXACT.u_s_mid{k}(end)^2 + ....
                     data_EXACT.v_s_mid{k}(end)^2))/2) * ...
                filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL));
% f5 = orrected omega
f5 = @(in_d, i)((1 - ((in_d.u_s_mid{i}(1) * in_d.x_s_mid{i}(2).^2 - ...
                       in_d.u_s_mid{i}(2) * in_d.x_s_mid{i}(1).^2) ./ ...
                      (in_d.x_s_mid{i}(2).^2 - in_d.x_s_mid{i}(1).^2) .* ...
                      in_d.F{i}).^2) * ...
                filter_trick((min(in_d.u_s_mid{i})) > MIN_HORZ_VEL));
fmass = @(in_d, i)(filter_trick(min(in_d.u_mid{i}) > MIN_HORZ_VEL) * ( ...
    ... % downstream component + upstream component
    trapz([-in_d.x_mid{i}(end:-1:1) in_d.x_mid{i}],  ...
          [in_d.eta_mid{i}(end:-1:1) in_d.eta_mid{i}])));

fcirc = @(in_d, i)(filter_trick(min(in_d.u_mid{i}) > MIN_HORZ_VEL) * (...
    ... % (phi + x)(infty) -- upstream component
    (-in_d.phi_mid{i}(end)+in_d.x_mid{i}(end)) - ...
    ... % (phi + x)(-infty) -- downstream component
    (in_d.phi_mid{i}(end)-in_d.x_mid{i}(end))));
fT = @(in_d, i)(0.5 * (fmass(in_d,i) - fcirc(in_d,i)) * ...
                  filter_trick(min(in_d.u_mid{i}) > MIN_HORZ_VEL));
fV = @(in_d, i)(filter_trick(min(in_d.u_mid{i}) > MIN_HORZ_VEL) * ( ...
     ... % downstream component, upstream component
     trapz([-in_d.x_mid{i}(end:-1:1), in_d.x_mid{i}], ...
           [in_d.eta_mid{i}(end:-1:1).*in_d.eta_mid{i}(end:-1:1), ...
            in_d.eta_mid{i}.*in_d.eta_mid{i}])/(2*in_d.F{i}^2)));
fE = @(in_d, i)(filter_trick(min(in_d.u_mid{i}) > MIN_HORZ_VEL) * ( ...
                  fT(in_d,i) + fV(in_d,i)) .* (in_d.F{i}^2));
                

j = 1;

c = blue_black(length(Avec));

for A = Avec
  old_wd = pwd();
  cd('../output/');
  try
    astruct = regexp(sprintf('%g', abs(A)), ...
                '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
    astr = strcat(fstruct.decimal);
    if length(astr) < 3
      astr = [astr, repmat('0', 1,3-length(astr))];
    end % zero padding
    
    evalstr = sprintf( ...
    ['if ~exist(''data_A%s'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''pressure_TYPE4_A%s_B280.m'', ''file'') || ...\n' ...
     '     exist(''pressure_TYPE4_A%s_B280.mat'', ''file'')\n' ...
     '    data_A%s = convert_and_load(''pressure_TYPE4_A%s_B280'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_A%s'', ''var'')\n' ...
     '  disp(''Calue of N : '');\n' ...
     '  disp(length(data_A%s.x_s_mid{1})+1);\n' ...
     '  data_A%s.branch = ''AH1'';\n' ...
     'end\n'], ...
     astr, astr, astr, astr, astr, astr, astr, astr);
    eval(evalstr)
    
  catch e
    cd(old_wd)
    throw(e)
  end
  
  cd(old_wd);
  
  this_color = c(j,:);
  
  figure(batch_fig_h)
  
  hold on
  
  evalstr = sprintf( ...
  ['if exist(''data_A%s'',''var'')\n' ...
   '  [~] = batch_plot_m2tikz(data_A%s, ...\n' ...
   '                            f5, f2,  ''-'', ...\n' ...
   '                            ''color'', this_color);\n' ...
... %   '[~, ~] = batch_plot(data_A%s, f5, f2, ''o'', this_color, 0.75);\n' ...
   'end'], ...
  astr, astr);
%   astr, astr, astr);

  eval(evalstr)
  
end



%% Choose Froude number and number of profiles
% Here a single Froude number is selected for further playing around. The
% plotting of the profiles and saving them to output is determined by a
% callback function (FSM).

disp(' ')
disp('---------------- Compile sets of profiles on a separate image(s)');
disp(' ')

INVALID_SELECTIONS = true;

while INVALID_SELECTIONS
  
  % Select Froude number
  disp_str1 = '';
  disp_str2 = '';

  for k = 1:length(Fvec)
    disp_str1 = [disp_str1, sprintf('%i\t', k)];
    disp_str2 = [disp_str2, sprintf('%.3f\t', Fvec(k))];
    if ~mod(k,10)
      disp(disp_str1);
      disp(disp_str2);
      disp(' ');
      disp_str1 = '';
      disp_str2 = '';
    end
  end
  disp(disp_str1);
  disp(disp_str2);
  
  prompt_str = input('Select Froude number (via index) : ','s');
  Findex = str2double(prompt_str);
  if ~ismember(Findex, 1:length(Fvec))
    disp('Invalid selection for F')
    continue
  end
  
  F = Fvec(Findex);
  
  fstruct = regexp(sprintf('%g', F), ...
           '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
  fstr = strcat(fstruct.digit, fstruct.decimal);
                     
  this_color = [0 0 0]; %c(1,:); 
  
  figure(branch_fig)
  
  hold on
  
  evalstr = sprintf( ...
  ['if exist(''data_F%s_STD'',''var'')\n' ...
   '  [~, ...\n' ...
   '   ~, ...\n' ...
   '   cb_xd, ...\n' ...
   '   cb_yd, ...\n' ...
   '   cbh] = batch_plot_m2tikz(data_F%s_STD, ...\n' ...
   '                                  f2, f3,  ''-'', ...\n' ...
   '                                  ''color'', this_color);\n' ...
   'for th=cbh\n' ...
   '  set(th,''UserData'', struct(''data'', data_F%s_STD, ... \n' ...
   '                              ''xd'', cb_xd, ''yd'', cb_yd));\n' ...
   'end\n' ...
   'end'], ...
   fstr, fstr, fstr);
 
  eval(evalstr)

  evalstr = sprintf( ...
  ['if exist(''data_F%s_DBL'',''var'')\n' ...
   '  [~, ...\n' ...
   '   ~, ...\n' ...
   '   cb_xd, ...\n' ...
   '   cb_yd, ...\n' ...
   '   cbh] = batch_plot_m2tikz(data_F%s_DBL, ...\n' ...
   '                                f2, f3,  ''-'', ...\n' ...
   '                                ''color'', this_color);\n' ...
   'for th=cbh\n' ...
   '  set(th,''UserData'', struct(''data'', data_F%s_DBL, ... \n' ...
   '                              ''xd'', cb_xd, ''yd'', cb_yd));\n' ...
   'end\n' ...
   'end'], ...
   fstr, fstr, fstr);
  eval(evalstr)
  
  evalstr = sprintf( ...
  ['if exist(''data_F%s_DIP'',''var'')\n' ...
   '  [~, ...\n' ...
   '   ~, ...\n' ...
   '   cb_xd, ...\n' ...
   '   cb_yd, ...\n' ...
   '   cbh] = batch_plot_m2tikz(data_F%s_DIP, ...\n' ...
   '                                f2, f3,  ''-'', ...\n' ...
   '                                ''color'', this_color);\n' ...
   'for th=cbh\n' ...
   '  set(th,''UserData'', struct(''data'', data_F%s_DIP, ... \n' ...
   '                              ''xd'', cb_xd, ''yd'', cb_yd));\n' ...
   'end\n' ...
   'end'], ...
   fstr, fstr, fstr);
 
  eval(evalstr)
  
  % This is where I should add in the u_min plot! This could take a while
  % though.
  
  % Ask user for values of U to use
  
  disp('----------------------- Display a branch with fixed u at crest');
  
  olddir = pwd();
  try
    cd '../output/'
    disp('Available fixed u branches (allegedly):');
    a = dir('pressure_U*_B280.m');
    b = dir('pressure_U*_B280.mat');
    a1 = cell(1,length(a));
    for k=1:length(a)
      a1{k} = getfield(a, {k}, 'name');
    end
    b1 = cell(1,length(b));
    for k=1:length(b)
      b1{k} = getfield(b, {k}, 'name');
    end

    a2 = regexp(a1,'_U(?<number>\d+)_', 'names');
    b2 = regexp(b1,'_U(?<number>\d+)_', 'names');

    a3 = zeros(1,length(a2));
    for k=1:length(a2)
      a3(k) = str2double(a2{k}.number);
    end
    b3 = zeros(1,length(b2));
    for k=1:length(b2)
      b3(k) = str2double(b2{k}.number);
    end

    U_avail = unique([a3, b3]);
    U_avail = sort(U_avail./100)

  catch err
    cd(olddir)
    throw(err)
  end
  cd(olddir)

  Uindex = input('Input u of interest BY INDEX (0 = none), Uindex = ');

  assert(isscalar(Uindex), 'invalid selection for U')
  assert(ismember(Uindex,0:length(U_avail)), 'invalid selection for u')

  if Uindex
    ustruct = regexp(sprintf('%0.2f', U_avail(Uindex)), ...
                     '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
    ustr = strcat('0', ustruct.decimal);
    old_wd = pwd();
    cd('../output/');
    try
      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_SL1'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''pressure_U%s_SL1_B280.m'', ''file'') || ...\n' ...
       '     exist(''pressure_U%s_SL1_B280.mat'', ''file'')\n' ...
       '    data_U%s_SL1 = convert_and_load(''pressure_U%s_SL1_B280'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_SL1'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, B = 2.80, SL1 soln'');\n' ...
       '  disp(length(data_U%s_SL1.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_SL2'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''pressure_U%s_SL2_B280.m'', ''file'') || ...\n' ...
       '     exist(''pressure_U%s_SL2_B280.mat'', ''file'')\n' ...
       '    data_U%s_SL2 = convert_and_load(''pressure_U%s_SL2_B280'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_SL2'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, B = 2.80, SL2 soln'');\n' ...
       '  disp(length(data_U%s_SL2.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_DL1'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''pressure_U%s_DL1_B280.m'', ''file'') || ...\n' ...
       '     exist(''pressure_U%s_DL1_B280.mat'', ''file'')\n' ...
       '    data_U%s_DL1 = convert_and_load(''pressure_U%s_DL1_B280'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_DL1'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, B = 2.80, DL1 soln'');\n' ...
       '  disp(length(data_U%s_DL1.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_DL2'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''pressure_U%s_DL2_B280.m'', ''file'') || ...\n' ...
       '     exist(''pressure_U%s_DL2_B280.mat'', ''file'')\n' ...
       '    data_U%s_DL2 = convert_and_load(''pressure_U%s_DL2_B280'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_DL2'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, B = 2.80, DL2 soln'');\n' ...
       '  disp(length(data_U%s_DL2.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

    catch err
      disp(err)
      cd(olddir)
      throw(err)
    end
    cd(olddir)

    evalstr = sprintf( ...
    ['if exist(''data_U%s_SL1'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_SL1, ...\n' ...
     '                        f2, f3,  ''--'', ...\n' ...
     '                        ''color'', [0 0 0]);\n' ...
     'end'], ...
     ustr, ustr);

    eval(evalstr)

    evalstr = sprintf( ...
    ['if exist(''data_U%s_SL2'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_SL2, ...\n' ...
     '                        f2, f3,  ''--'', ...\n' ...
     '                        ''color'', [0 0 0]);\n' ...
     'end'], ...
     ustr, ustr);

    eval(evalstr)

    evalstr = sprintf( ...
    ['if exist(''data_U%s_DL1'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_DL1, ...\n' ...
     '                        f2, f3,  ''--'', ...\n' ...
     '                        ''color'', [1 0 0]);\n' ...
     'end'], ...
     ustr, ustr);

    eval(evalstr)

    evalstr = sprintf( ...
    ['if exist(''data_U%s_DL2'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_DL2, ...\n' ...
     '                        f2, f3,  ''--'', ...\n' ...
     '                        ''color'', [0 0 1]);\n' ...
     'end'], ...
     ustr, ustr);

    eval(evalstr)
  else
    disp('Not displaying any fixed u branch')
  end

  num_profiles = 0/0;
  while isnan(num_profiles)
    prompt_str = input('Total # profiles to select (1) : ', 's');
    if strcmp(prompt_str, '')
      num_profiles = 1;
    else
      num_profiles = str2double(prompt_str);
    end
  end

  % Display selected parameters
  disp(' ')
  disp('----------------------- Summary of parameters chosen for overlay'); 
  disp(['F : ' fstr ]);
  disp(['total #profiles : ' num2str(num_profiles)]);
  disp(' ')

  % Check that this is all ok!
  prompt_str = input('Proceed (y)/n : ', 's');
  if ~strcmp(prompt_str, 'n')
    INVALID_SELECTIONS = false;
  end

end


set(batch_fig_h, 'BusyAction', 'cancel', 'Interruptible', 'off')
  
evalstr = sprintf(['set(batch_fig_h, ''WindowButtonUpFcn'', ...\n' ...
                   '     {@callback_pressure_fixA_m2tikz, ...\n' ...
                   '                                   F, ...\n' ...
                   '                         profile_fig, ...\n' ...
                   '                           theta_fig, ...\n' ...
                   '                        branch_fig_h, ...\n' ...
                   '                        branch_fig_E, ...\n' ...
                   '                        num_profiles});\n']);
            
eval(evalstr)

disp('Ready to select branch (click figure)');
      
%figure(batch_fig_h)

% Need to add option to do custom axis!
% sigh ...
% And add code that changes the ticklabels to the right font!
%axis tight

%bf_axis = axis();

%bf_axis(1) = 0.992; %0.85;
%bf_axis(2) = 1.0; %
%bf_axis(3) = 1.29085;
%bf_axis(4) = 1.29090;%
%axis(bf_axis)

%stretch(gcf,[1.05 1.1])

%bf_axis = axis();


%ylabel('F', 'rotation', 0, 'UserData','matlabfrag:$F$');
%xlabel('w', 'UserData','matlabfrag:$\omega$');
%ylabel('$F$')
%xlabel('$\omega$')
%grid on
%box on

% Need to enter a name for this mother

%if PRINT_GRAPHS_TO_FILE
%  cd('../output/')
%  file_name = strcat(date_str, ...
%                     '_pressure_B280_STO', ...
%                     '.tikz');
%  matlab2tikz(file_name, 'width','\figurewidth','parseStrings',false)
 % matlabfrag(file_name, 'handle', batch_fig_h);
%  cd('../matlab/')
%end


%figure(branch_fig)

% Fstr = num2str(round(F*100));
%axis(bf_axis)

%ylabel('F', 'rotation', 0, 'UserData','matlabfrag:$F$');
%xlabel('w', 'UserData','matlabfrag:$\omega$');
%ylabel('$F$')
%xlabel('$\omega$')

%grid on
%box on

% Need to enter a name for this mother

%if PRINT_GRAPHS_TO_FILE
%  cd('../output/')
%  Astr = sprintf('%03i', round(A*100));

%  file_name = sprintf('%s_pressure_A%s_B280_height_summary.tikz', ...
%                      date_str, Astr);
%  matlab2tikz(file_name, 'width','\figurewidth','parseStrings',false)
  %matlabfrag(file_name, 'handle', branch_fig);
%  cd('../matlab/')

%end

