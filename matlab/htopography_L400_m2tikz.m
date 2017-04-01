%% Create plots of $$ A-\eta(0) $$ plane and profile subfigures
%
% To-do: Add in plot of fixed u_min.
%
% Stephen Wade 12/12/2013

global PRINT_GRAPHS_TO_FILE;
global AXES_LIMITS;

AXES_LIMITS = [-0.22, 0.22, -0.1, 1.0];

% Changed naming convention for files to improve sorting.
BATCHPLOT_STRETCH = 1.15;
% self explanatory
PRINT_GRAPHS_TO_FILE = true;
% This is the velocity at which 'stagnation' is assumed to occur
MIN_HORZ_VEL = 0.025;

disp('------------------------------ Process Smoothed-box data : L = 4.00');
disp(' ');

prompt_str = input('Recompute results y/(n)? : ','s');
if strcmp(prompt_str, 'n') || strcmp(prompt_str, '')
  RECOMPUTE = false;
else
  RECOMPUTE = true;
end

if RECOMPUTE
  clear -regexp cb_xd|cb_yd|data_F|data_U|xd_|yd_
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
  disp('Available Froude numbers (allegedly):');
  a = dir('htopography_F*_L400.m');
  b = dir('htopography_F*_L400.mat');
  a1 = cell(1,length(a));
  for k=1:length(a)
    a1{k} = getfield(a, {k}, 'name');
  end
  b1 = cell(1,length(b));
  for k=1:length(b)
    b1{k} = getfield(b, {k}, 'name');
  end
  
  a2 = regexp(a1,'_F(?<number>\d+)_', 'names');
  b2 = regexp(b1,'_F(?<number>\d+)_', 'names');
  
  a3 = zeros(1,length(a2));
  for k=1:length(a2)
    a3(k) = str2double(a2{k}.number);
  end
  b3 = zeros(1,length(b2));
  for k=1:length(b2)
    b3(k) = str2double(b2{k}.number);
  end
  
  F_avail = unique([a3, b3]);
  F_avail = sort(F_avail.*(10.^(-floor(log10(F_avail)))))
  
catch err
  cd(olddir)
  throw(err)
end
cd(olddir)

Findex = input('Input Froude #s of interest BY INDEX, Findex = ');

assert(isvector(Findex), 'invalid selection for F')
assert(sum(ismember(Findex,1:length(F_avail))) == length(Findex), ...
       'invalid selection for F')

Fvec = F_avail(Findex);

set(0,'DefaultTextInterpreter','none');


%% Let's get to it.

disp('Ok, valid selection, let''s get to it')

batch_fig_h    = figure('units', 'centimeters', ...
                        'position', [0.5 9 7 5.25]);
                      
profile_fig = figure('units', 'centimeters', ...
                     'position', [10.5 9 7 5.25]);

branch_fig = figure('units', 'centimeters', ...
                    'position', [20.5 9 7 5.25]);

% Possible filter tricks
%
% 1. Check the minimum velocity (for stagnation)
% filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL)

% use omega
%f1 = @(in_d, i)((1 - (in_d.F{i} * min(in_d.u_s_mid{i}))^2) * ...
%                filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL));

% use A
%f2 = @(in_d, i)((in_d.A{i}) *  ...
%               filter_trick(exp(min(in_d.tau_s_mid{i})) > ...
%                            MIN_HORZ_VEL));

% use AREA
f2 = @(in_d, i)((2*trapz(in_d.x_b_mid{i}, in_d.y_b_mid{i})) *  ...
               filter_trick(exp(min(in_d.tau_s_mid{i})) > ...
                            MIN_HORZ_VEL));


% use eta0
f3 = @(in_d, i)((in_d.eta_s_mid{i}(1) + ...
                 ((in_d.F{i}^2) * ...
                  (1 - exp(in_d.tau_s_mid{i}(end)) .* ...
                   cos(in_d.theta_s_mid{i}(end))))) * ...
                filter_trick(exp(min(in_d.tau_s_mid{i})) > ...
                             MIN_HORZ_VEL));
% use F
%f4 = @(in_d, i)(in_d.F{i} * ...
%                filter_trick(min(in_d.u_s_mid{i}) > MIN_HORZ_VEL));

j = 1;

c = blue_black(length(Fvec));

for F = Fvec
  old_wd = pwd();
  cd('../output/');
  try
    fstruct = regexp(sprintf('%g', F), ...
                '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
    if length(fstruct.decimal) < 2
      fstruct.decimal = strcat(fstruct.decimal, ...
                               repmat('0', 1, 2-length(fstruct.decimal)));
    end
    fstr = strcat(fstruct.digit, fstruct.decimal);

    % Load data using some annoying code.
    
    evalstr = sprintf( ...
    ['if ~exist(''data_F%s_DBL'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''htopography_F%s_DBL_L400.m'', ''file'') || ...\n' ...
     '     exist(''htopography_F%s_DBL_L400.mat'', ''file'')\n' ...
     '    data_F%s_DBL = convert_and_load(''htopography_F%s_DBL_L400'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_F%s_DBL'', ''var'')\n' ...
     '  disp(''Value of N, for F=%f, L = 4.00, DBL soln'');\n' ...
     '  disp(length(data_F%s_DBL.x_s_mid{1})+1);\n' ...
     '  data_F%s_DBL.branch = ''DBL'';\n' ...
     'end\n'], ...
     fstr, fstr, fstr, fstr, fstr, fstr, F, fstr, fstr);
    eval(evalstr)

    evalstr = sprintf( ...
    ['if ~exist(''data_F%s_DIP'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''htopography_F%s_US2_L400.m'', ''file'') || ...\n' ...
     '     exist(''htopography_F%s_US2_L400.mat'', ''file'')\n' ...
     '    data_F%s_DIP = convert_and_load(''htopography_F%s_US2_L400'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_F%s_DIP'', ''var'')\n' ...
     '  disp(''Value of N, for F=%f, L = 4.00, DIP soln'');\n' ...
     '  disp(length(data_F%s_DIP.x_s_mid{1})+1);\n' ...
     '  data_F%s_DIP.branch = ''DIP'';\n' ...
     'end\n'], ...
     fstr, fstr, fstr, fstr, fstr, fstr, F, fstr, fstr);
   
    eval(evalstr)
   
    evalstr = sprintf( ...
    ['if ~exist(''data_F%s_STD'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''htopography_F%s_US1_L400.m'', ''file'') || ...\n' ...
     '     exist(''htopography_F%s_US1_L400.mat'', ''file'')\n' ...
     '    data_F%s_STD = convert_and_load(''htopography_F%s_US1_L400'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_F%s_STD'', ''var'')\n' ...
     '  disp(''Value of N, for F=%f, L = 4.00, STD soln'');\n' ...
     '  disp(length(data_F%s_STD.x_s_mid{1})+1);\n' ...
     '  data_F%s_STD.branch = ''STD'';\n' ...
     'end\n'], ...
     fstr, fstr, fstr, fstr, fstr, fstr, F, fstr, fstr);
   
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
  ['if exist(''data_F%s_STD'',''var'')\n' ...
   '  [~] = batch_plot_m2tikz(data_F%s_STD, ...\n' ...
   '                            f2, f3,  ''-'', ...\n' ...
   '                            ''color'', this_color);\n' ...
... %   '[~, ~] = batch_plot(data_F%s_STD, f2, f3, ''o'', this_color, 0.75);\n' ...
   'end'], ...
  fstr, fstr);
%   fstr, fstr, fstr);

  eval(evalstr)

  evalstr = sprintf( ...
  ['if exist(''data_F%s_DBL'',''var'')\n' ...
   '  [~] = batch_plot_m2tikz(data_F%s_DBL, ...\n' ...
   '                          f2, f3,  ''-'', ...\n' ...
   '                          ''color'', this_color);\n' ...
... %   '[~, ~] = batch_plot(data_F%s_DBL, f2, f3, ''o'', this_color, 0.75);\n' ...
   'end'], ...
   fstr, fstr);
%   fstr, fstr, fstr);
 
  eval(evalstr)
  
  evalstr = sprintf( ...
  ['if exist(''data_F%s_DIP'',''var'')\n' ...
   '[~] = batch_plot_m2tikz(data_F%s_DIP, ...\n' ...
   '                        f2, f3,  ''-'', ...\n' ...
   '                        ''color'', this_color);\n' ...
... %   '[~, ~] = batch_plot(data_F%s_DIP, f2, f3, ''o'', this_color, 0.75);\n' ...
   'end'], ...
   fstr, fstr);
%   fstr, fstr, fstr);
 
  eval(evalstr)
  
  j = j + 1;
  
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
  if length(fstruct.decimal) < 2
    fstruct.decimal = strcat(fstruct.decimal, ...
                             repmat('0', 1, 2-length(fstruct.decimal)));
  end
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
  
  disp('------------------------- Display a branch with fixed u at crest');
  
  olddir = pwd();
  try
    cd '../output/'
    disp('Available fixed u branches (allegedly):');
    a = dir('htopography_U*_L400.m');
    b = dir('htopography_U*_L400.mat');
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
    if isempty(U_avail)
      disp('  None.')
    else
      U_avail = sort(U_avail./100)
    end
    
  catch err
    cd(olddir)
    throw(err)
  end
  cd(olddir)
  
  if ~isempty(U_avail)

    Uindex = input('Input u of interest BY INDEX, Uindex = ');

    assert(isscalar(Uindex), 'invalid selection for U')
    assert(ismember(Uindex,1:length(U_avail)), 'invalid selection for u')

    ustruct = regexp(sprintf('%0.2f', U_avail(Uindex)), ...
                     '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
    ustr = strcat('0', ustruct.decimal);
    old_wd = pwd();
    cd('../output/');
    try
      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_SS1'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''htopography_U%s_SS1_L400.m'', ''file'') || ...\n' ...
       '     exist(''htopography_U%s_SS1_L400.mat'', ''file'')\n' ...
       '    data_U%s_SS1 = convert_and_load(''htopography_U%s_SS1_L400'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_SS1'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, L = 4.00, SS1 soln'');\n' ...
       '  disp(length(data_U%s_SS1.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_SS2'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''htopography_U%s_SS2_L400.m'', ''file'') || ...\n' ...
       '     exist(''htopography_U%s_SS2_L400.mat'', ''file'')\n' ...
       '    data_U%s_SS2 = convert_and_load(''htopography_U%s_SS2_L400'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_SS2'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, L = 4.00, SS2 soln'');\n' ...
       '  disp(length(data_U%s_SS2.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)

      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_DL1'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''htopography_U%s_DL1_L400.m'', ''file'') || ...\n' ...
       '     exist(''htopography_U%s_DL1_L400.mat'', ''file'')\n' ...
       '    data_U%s_DL1 = convert_and_load(''htopography_U%s_DL1_L400'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_DL1'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, L = 4.00, DL1 soln'');\n' ...
       '  disp(length(data_U%s_DL1.x_s_mid{1})+1);\n' ...
       'end\n'], ...
       ustr, ustr, ustr, ustr, ustr, ustr, U_avail(Uindex), ustr);
      eval(evalstr)
      
      evalstr = sprintf( ...
      ['if ~exist(''data_U%s_DL2'', ''var'') || RECOMPUTE\n' ...
       '  if exist(''htopography_U%s_DL2_L400.m'', ''file'') || ...\n' ...
       '     exist(''htopography_U%s_DL2_L400.mat'', ''file'')\n' ...
       '    data_U%s_DL2 = convert_and_load(''htopography_U%s_DL2_L400'');\n' ...
       '  end\n' ...
       'end\n' ...
       'if exist(''data_U%s_DL2'', ''var'')\n' ...
       '  disp(''Value of N, for U=%f, L = 4.00, DL2 soln'');\n' ...
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
    ['if exist(''data_U%s_SS1'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_SS1, ...\n' ...
     '                        f2, f3,  ''--'', ...\n' ...
     '                        ''color'', [0 0 0]);\n' ...
     'end'], ...
     ustr, ustr);

    eval(evalstr)

    evalstr = sprintf( ...
    ['if exist(''data_U%s_SS2'',''var'')\n' ...
     '[~] = batch_plot_m2tikz(data_U%s_SS2, ...\n' ...
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
                   '         {@callback_htopography_m2tikz, ...\n' ...
                   '                                     F, ...\n' ...
                   '                           profile_fig, ...\n' ...
                   '                            branch_fig, ...\n' ...
                   '                          num_profiles});\n']);
            
eval(evalstr)

disp('Ready to select branch (click figure)');

