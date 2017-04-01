%% Print out the type 3 solutions to hyperbolic function topography
%
% Stephen Wade 15/07/2016
%
% Adapted from htopography_TYPE2_m2tikz
%
% Choose a value of A (which will be fixed) and print out the energy and
% other various thingos as omega increases towards one.

global PRINT_GRAPHS_TO_FILE;
global AXES_LIMITS;
global AXES_LIMITS_MASS;

AXES_LIMITS      = [0.7, 1.00, 1.25, 1.34]; % will just modify these in tikz
AXES_LIMITS_MASS = [0.7, 1.00, 0, 2.0];   % will just modify these in tikz

% Changed naming convention for files to improve sorting.
BATCHPLOT_STRETCH = 1.15;
% self explanatoryy
PRINT_GRAPHS_TO_FILE = true;
% This is the velocity at which 'stagnation' is assumed to occur
MIN_HORZ_VEL = 0.02; % works ok enough for F=1.262, 0.02 for others
% note: must be aware of make change to batch_plot_m2tikz for F=1.262 plot
WAVE_TYPE = 3;

ustruct = regexp(sprintf('%g', MIN_HORZ_VEL), ...
                 '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
ustr = strcat(ustruct.digit, ustruct.decimal);
if length(ustr) < 3
  ustr = [ustr, repmat('0', 1, 3-length(ustr))];
end % zero padding

disp('------------------ Process hyperbolic topography : TYPE 3 Solution');
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

if exist('branch_fig_mass','var')
  if ishandle(branch_fig_mass)
    close(branch_fig_mass)
  end
  clear branch_fig_mass
end

if exist('branch_fig_h','var')
  if ishandle(branch_fig_h)
    close(branch_fig_h)
  end
  clear branch_fig_h
end

date_str = strrep(date(),'-','_');

% Ask user for values of amplitude of forcing, A, to use
olddir = pwd();
try
  cd '../output/'
  disp('Available amplitude of forcing (allegedly):');
  a = dir(['htopography_A*_U', ustr, '_TYPE3.m']);
  b = dir(['htopography_A*_U', ustr, '_TYPE3.mat']);
  a1 = cell(1,length(a));
  for k=1:length(a)
    a1{k} = getfield(a, {k}, 'name');
  end
  b1 = cell(1,length(b));
  for k=1:length(b)
    b1{k} = getfield(b, {k}, 'name');
  end
  
  a2 = regexp(a1,'_A(?<number>\d+)_', 'names');
  b2 = regexp(b1,'_A(?<number>\d+)_', 'names');
  
  a3 = zeros(1,length(a2));
  for k=1:length(a2)
    a3(k) = str2double(a2{k}.number);
  end
  b3 = zeros(1,length(b2));
  for k=1:length(b2)
    b3(k) = str2double(b2{k}.number);
  end

  A_avail = unique([a3, b3]);
  A_avail = sort(A_avail.*(10.^(-2)))
  
catch err
  cd(olddir)
  throw(err)
end
cd(olddir)

Aindex = input('Input A of interest BY INDEX : ');

assert(isvector(Aindex), 'invalid selection for A')
assert(sum(ismember(Aindex,1:length(A_avail))) == length(Aindex), ...
       'invalid selection for A');

Avec = A_avail(Aindex);

set(0,'DefaultTextInterpreter','none');

%% Let's get to it.

disp('Ok, valid selection, let''s get to it')

batch_fig_h  = figure('units', 'centimeters', ...
                       'position', [0.5 9 7 5.25]);
profile_fig  = figure('units', 'centimeters', ...
                      'position', [10.5 9 7 5.25]);
theta_fig    = figure('units', 'centimeters', ...
                      'position', [20.5 9 7 5.25]);
branch_fig_h    = figure('units', 'centimeters', ...
                         'position', [30.5 9 7 5.25]);
branch_fig_mass = figure('units', 'centimeters', ...
                         'position', [40.5 9 7 5.25]);

% f1 = omega
f1 = @(in_d, i)((1 - (in_d.FROUDE{i} * min(exp(in_d.tau_s_mid{i})))^2) * ...
                filter_trick((min(exp(in_d.tau_s_mid{i}))) > MIN_HORZ_VEL));
% f2 = F
f2 = @(in_d, i)((in_d.FROUDE{i}) *  ...
                filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL));
% f3 = corrected eta0
f3 = @(in_d, i)((in_d.eta_s_mid{i}(1) + ...
                 (in_d.FROUDE{i}^2) * ...
                 (1-(exp(in_d.tau_s_mid{i}(end))^2 + ....
                     in_d.v_s_mid{i}(end)^2))/2) * ...
                filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL));
% f5 = corrected omega - for crests at 0!
f5 = @(in_d, i)((1 - ((exp(in_d.tau_s_mid{i}(1)) * in_d.x_s_mid{i}(2) - ...
                       exp(in_d.tau_s_mid{i}(2)) * in_d.x_s_mid{i}(1)) ./ ...
                      (in_d.x_s_mid{i}(2) - in_d.x_s_mid{i}(1)) .* ...
                      in_d.FROUDE{i}).^2) * ...
                filter_trick((min(exp(in_d.tau_s_mid{i}))) > MIN_HORZ_VEL));
fmass = @(in_d, i)(filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL) * ( ...
    ... % downstream component + upstream component
    trapz([-in_d.x_s_mid{i}(end:-1:1) in_d.x_s_mid{i}],  ...
          [in_d.eta_s_mid{i}(end:-1:1) in_d.eta_s_mid{i}])));

fcirc = @(in_d, i)(filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL) * (...
    ... % (phi + x)(infty) -- upstream component
    (-in_d.phi_s_mid{i}(end)+in_d.x_s_mid{i}(end)) - ...
    ... % (phi + x)(-infty) -- downstream component
    (in_d.phi_s_mid{i}(end)-in_d.x_s_mid{i}(end))));
fT = @(in_d, i)(0.5 * (fmass(in_d,i) - fcirc(in_d,i)) * ...
                  filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL));
fV = @(in_d, i)(filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL) * ( ...
     ... % downstream component, upstream component
     trapz([-in_d.x_s_mid{i}(end:-1:1), in_d.x_s_mid{i}], ...
           [in_d.eta_s_mid{i}(end:-1:1).*in_d.eta_s_mid{i}(end:-1:1), ...
            in_d.eta_s_mid{i}.*in_d.eta_s_mid{i}])/(2*in_d.FROUDE{i}^2)));
fE = @(in_d, i)(filter_trick(min(exp(in_d.tau_s_mid{i})) > MIN_HORZ_VEL) * ( ...
                  fT(in_d,i) + fV(in_d,i)) .* (in_d.FROUDE{i}^2));
                

j = 1;

c = blue_black(length(Avec));

for A = Avec
  old_wd = pwd();
  cd('../output/');
  try
    astruct = regexp(sprintf('%g', A), ...
                     '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
    astr = strcat(astruct.digit,astruct.decimal);
    if length(astr) < 3
      astr = [astr, repmat('0', 1, 3-length(astr))];
    end % zero padding
    temp_list = {astr, ...
                 astr, ustr, ...
                 astr, ustr, ...
                 astr, astr, ustr, ...
                 astr, ...
                 astr, ...
                 astr};
    evalstr = sprintf( ...
    ['if ~exist(''data_A%s'', ''var'') || RECOMPUTE\n' ...
     '  if exist(''htopography_A%s_U%s_TYPE3.m'', ''file'') || ...\n' ...
     '     exist(''htopography_A%s_U%s_TYPE3.mat'', ''file'')\n' ...
     '    data_A%s = convert_and_load(''htopography_A%s_U%s_TYPE3'');\n' ...
     '  end\n' ...
     'end\n' ...
     'if exist(''data_A%s'', ''var'')\n' ...
     '  disp(''Calue of N : '');\n' ...
     '  disp(length(data_A%s.x_s_mid{1})+1);\n' ...
     '  data_A%s.branch = ''AH1'';\n' ...
     'end\n'], ...
     temp_list{:});
     eval(evalstr)
    
  catch e
    cd(old_wd)
    throw(e)
  end
  
  cd(old_wd);
  
  this_color = c(j,:);
  
  figure(batch_fig_h)
  
  hold on
%  temp_list = {astr, astr};
  temp_list = {astr, astr, astr};
  evalstr = sprintf( ...
  ['if exist(''data_A%s'',''var'')\n' ...
   '  [~] = batch_plot_m2tikz(data_A%s, ...\n' ...
   '                          f1, f2,  ''-'', ...\n' ...
   '                          ''color'', this_color);\n' ...
   '[~, ~] = batch_plot(data_A%s, f1, f2, ''o'', this_color, 0.75);\n' ...
   'end'], ...
   temp_list{:});
 
  eval(evalstr)
  
end

%% Choose A and number of profiles
% Here a single A is selected for further playing around. The
% plotting of the profiles and saving them to output is determined by a
% callback function (FSM).
% 29/07/2014             - what the hell does FSM mean? Flying Spaghetti Monster?
% 29/07/2014 + 5 seconds - oh, it means Finite State Machine. 

disp(' ')
disp('---------------- Compile sets of profiles on a separate image(s)');
disp(' ')

INVALID_SELECTIONS = true;

while INVALID_SELECTIONS
  
  % Select amplitude of forcing
  disp_str1 = '';
  disp_str2 = '';

  for k = 1:length(Avec)
    disp_str1 = [disp_str1, sprintf('%i\t', k)];
    disp_str2 = [disp_str2, sprintf('%.3f\t', Avec(k))];
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
  
  prompt_str = input('Select A (via index) : ','s');
  Aindex = str2double(prompt_str);
  if ~ismember(Aindex, 1:length(Avec))
    disp('Invalid selection for A')
    continue
  end
  
  A = Avec(Aindex);
  
  astruct = regexp(sprintf('%g', A), ...
                   '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
  astr = strcat(astruct.digit,astruct.decimal);
  if length(astr) < 3
    astr = [astr, repmat('0', 1, 3-length(astr))];
  end % zero padding
                     
  this_color = [0 0 0]; %c(1,:); 
  
  figure(branch_fig_h)
  
  hold on
  
  temp_list = {astr, astr, astr};
  
  evalstr = sprintf( ...
  ['if exist(''data_A%s'',''var'')\n' ...
   '  [~, ...\n' ...
   '   ~, ...\n' ...
   '   cb_xd, ...\n' ...
   '   cb_yd, ...\n' ...
   '   cbh] = batch_plot_m2tikz(data_A%s, ...\n' ...
   '                            f1, f2,  ''-'', ...\n' ...
   '                            ''color'', this_color);\n' ...
   'for th=cbh\n' ...
   '  set(th,''UserData'', struct(''data'', data_A%s, ... \n' ...
   '                              ''xd'', cb_xd, ''yd'', cb_yd));\n' ...
   'end\n' ...
   'end'], ...
   temp_list{:});
 
  eval(evalstr)
 
  hold off
 
  figure(branch_fig_mass)
  
  hold on
  
  evalstr = sprintf( ...
  ['if exist(''data_A%s'',''var'')\n' ...
   '  [~, ...\n' ...
   '   ~, ...\n' ...
   '   cb_xd, ...\n' ...
   '   cb_yd, ...\n' ...
   '   cbh] = batch_plot_m2tikz(data_A%s, ...\n' ...
   '                            f1, fmass,  ''-'', ...\n' ...
   '                            ''color'', this_color);\n' ...
   'for th=cbh\n' ...
   '  set(th,''UserData'', struct(''data'', data_A%s, ... \n' ...
   '                              ''xd'', cb_xd, ''yd'', cb_yd));\n' ...
   'end\n' ...
   'end'], ...
   temp_list{:});
 
  eval(evalstr)
  
  hold off
 
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
  disp(['A : ' astr ]);
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
                   '     {@callback_htopography_fixA_m2tikz, ...\n' ...
                   '                                      A, ...\n' ...
                   '                              WAVE_TYPE, ...\n' ...
                   '                            profile_fig, ...\n' ...
                   '                              theta_fig, ...\n' ...
                   '                           branch_fig_h, ...\n' ...
                   '                        branch_fig_mass, ...\n' ...
                   '                           num_profiles});\n']);
            
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

