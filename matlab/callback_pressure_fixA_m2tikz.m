function callback_pressure_fixA_m2tikz(src, ...
                                       ~, ...
                                       A, ...
                                       prof_fh, ...
                                       theta_fh, ...
                                       branch_h_fh, ...
                                       branch_E_fh, ...
                                       total_profiles)
% Callback for examining Gaussian pressure results.
% Using matlab2tikz
%
% todo : implement a choice about profile picture axis.
%
% written by Stephen Wade
% 17/12/2013

% 12/12/2013 Adapted from callback_overlay_pressure
% 21/5/2013  Original code adapted from callback_overlay_gtopography.m

  persistent CB_M2TIKZ_STATE
  persistent SOLN_BRANCH
  persistent PLOT_TOPOGRAPHY
  persistent PROF_PLOTFILE_NAME
  persistent THETA_PLOTFILE_NAME
  persistent SUB_PROFILE_COUNT
  persistent ALL_PROFILE_COUNT
  persistent NSUB_PROFILES
  
  global AXES_LIMITS
  global AXES_LIMITS_E
  global PRINT_GRAPHS_TO_FILE
  
  DATE_STR = strrep(date(),'-','_');
  PROFPLOT_STRETCH = 1.1;
  
  % need some choice about profile picture axis limits.

  if isempty(CB_M2TIKZ_STATE)
    CB_M2TIKZ_STATE = 'CHOOSE_BRANCH';
  end
  
  if isempty(ALL_PROFILE_COUNT)
    ALL_PROFILE_COUNT = 0;
  end
  
  DEFAULT_PLOTSTYLE = '-';
  
  astruct = regexp(sprintf('%g', abs(A)), ...
                   '(?<digit>\d+)\.(?<decimal>\d+)', 'names');
  astr = strcat(astruct.decimal);
  if length(astr) < 3
    astr = [astr, repmat('0', 1,3-length(astr))];
  end
    
  switch CB_M2TIKZ_STATE
  case 'CHOOSE_BRANCH'
    % figure out what branch the user wants
    
    % first, get the data
    branch_axh = get(branch_h_fh, 'Children');
    branch_ax = get(branch_axh);
    ud = cell([1,length(branch_ax.Children)]);
    k = 1;
    for branch_ph = branch_ax.Children'
      ud{k} = get(branch_ph, 'UserData');
      k = k + 1;
    end
    
    disp(' ')
    disp('--------------------------------- Plot a selection of profiles');
    % Select dip loop, middle loop, or double loop?

    bcount = 0;
    for k = 1:length(ud)
      if ~isempty(ud{k})
        bcount = bcount+1;
      end
    end
    
    s = cell([1,bcount]);
    j = 1;
    for k = 1:length(ud)
      if ~isempty(ud{k})
        s{j} = ud{k}.data.branch;
        j = j + 1;
      end
    end
    
    s = unique(s);
    
    disp('------------------------------------------- Available branches');
    
   
    disp(s')

    prompt_str = input(['Select branch (' s{1} ') : '],'s');
    if ismember(prompt_str, s)
      SOLN_BRANCH = prompt_str;
    elseif strcmp(prompt_str, '')
      SOLN_BRANCH = s{1};
    else
      SOLN_BRANCH = [];
      error('Invalid branch selection');
    end
    
    prompt_str = input('Print topography+interior y/(n) : ','s');
    if ismember(lower(prompt_str), {'y', 'n'})
      if lower(prompt_str) == 'y'
        PLOT_TOPOGRAPHY = true;
      else
        PLOT_TOPOGRAPHY = false;
      end
    elseif strcmp(prompt_str, '')
        PLOT_TOPOGRAPHY = false;
    else
      PLOT_TOPOGRAPHY = [];
      error('Invalid selection');
    end

    % Choose number of profiles to select from this branch.
    prompt_str = input('Input # profiles to select from branch (1) : ', 's');
    if strcmp(prompt_str, '')
      NSUB_PROFILES = 1;
    else
      NSUB_PROFILES = str2double(prompt_str);
    end
    
    assert(~isnan(NSUB_PROFILES), 'Invalid input for # profiles');

    % Lets get a name for the plot of this stage.
    additional_str = sprintf('A%s_B280_xxxx_prof', astr);
    prompt_str = input(['Plot file name ''extension'' (', ...
                         additional_str, ') : ' ...
                         DATE_STR, '_pressure_'], 's');
    if ~strcmp(prompt_str, '')
      additional_str = prompt_str;
    end

    PROF_PLOTFILE_NAME = strcat(DATE_STR, ...
                                '_pressure_', ...
                                additional_str);

    % Lets get a name for the plot of this stage.
    additional_str = sprintf('A%s_B280_xxxx_theta', astr);
    prompt_str = input(['Plot file name ''extension'' (', ...
                         additional_str, ') : ' ...
                         DATE_STR, '_pressure_'], 's');
    if ~strcmp(prompt_str, '')
      additional_str = prompt_str;
    end

    THETA_PLOTFILE_NAME = strcat(DATE_STR, ...
                                 '_pressure_', ...
                                 additional_str);
    % Display selected parameters.
    disp(' ')
    disp('--------------------------------- Summary of parameters chosen'); 
    disp(['branch          : ' SOLN_BRANCH]);
    disp(['#profiles       : ' num2str(NSUB_PROFILES) ]);
    disp(['prof file name  : ' PROF_PLOTFILE_NAME ]);
    disp(['theta file name : ' THETA_PLOTFILE_NAME ]);


    % Check that this is all ok!
    prompt_str = input('Proceed (y)/n : ', 's');
    if ~strcmp(prompt_str, 'y') && ~strcmp(prompt_str, '')
      error('Invalid choice of parameters for profile selection')
    end

    % Option to clear?
    prompt_str = input('Clear profile pic? (y)/n : ', 's');
    if strcmp(prompt_str, 'y') || strcmp(prompt_str, '')
      clf(prof_fh)
    end
    
    SUB_PROFILE_COUNT = 0;
    CB_M2TIKZ_STATE = 'CHOOSE_PROFILE';
      
  case 'CHOOSE_PROFILE'
    % user has selected a profile, plot it
    
    % get the branch data
    branch_axh = get(branch_h_fh, 'Children');
    branch_ax = get(branch_axh);

    branch_E_axh = get(branch_E_fh, 'Children');
    branch_E_ax = get(branch_E_axh);
    
    FOUND_BRANCH = false;
    for branch_ph = branch_ax.Children'
      in_ud = get(branch_ph, 'UserData');
      if ~isempty(in_ud)
        in_d = in_ud.data;
        if isfield(in_d, 'branch')
          if strcmp(in_d.branch, SOLN_BRANCH)
            xd = in_ud.xd;
            yd = in_ud.yd;
            FOUND_BRANCH = true;
            break
          end
        end
      end
    end
    
    FOUND_E_BRANCH = false;
    for branch_E_ph = branch_E_ax.Children'
      in_E_ud = get(branch_E_ph, 'UserData');
      if ~isempty(in_E_ud)
        in_E_d = in_E_ud.data;
        if isfield(in_E_d, 'branch')
          if strcmp(in_E_d.branch, SOLN_BRANCH)
            xd_E = in_E_ud.xd;
            yd_E = in_E_ud.yd;
            FOUND_E_BRANCH = true;
            break
          end
        end
      end
    end
    
    if ~(FOUND_BRANCH && FOUND_E_BRANCH)
      error('Couldn''t find branch?')
    end
    
    cursor_position = get(src, 'CurrentPoint');
    fig_position = get(src, 'Position');
    axes_position = get(get(src, 'CurrentAxes'), 'Position');

    fig_width = fig_position(3);
    fig_height = fig_position(4);

    axes_width = fig_width * axes_position(3);
    axes_height = fig_height * axes_position(4);
  
    cursor_x = cursor_position(1) - (fig_width * axes_position(1));
    cursor_y = cursor_position(2) - (fig_height * axes_position(2));

    if cursor_x < 0 || cursor_x > axes_width || ...
       cursor_y < 0 || cursor_y > axes_height
      error('Can''t find nearby data')
    end

    c_axis = axis();

    xd_axes = axes_width * (xd - c_axis(1)) / (c_axis(2) - c_axis(1));
    yd_axes = axes_height * (yd - c_axis(3)) /(c_axis(4) - c_axis(3));

    [waste, i] = min((xd_axes - cursor_x).^2 + (yd_axes - cursor_y).^2);
    if waste < (axes_width*axes_height*0.01)

    disp('---------------------------- Overlaying gaussian pressure data');
    disp(' ');
    disp(['Selected profile number ' num2str(i)]);
    
    % need to add a 'colour/number' option
    
    disp(' ');
    disp('------------------------------ Parameters for graphical output');
    disp(' ');
    
    
    DEFAULT_COLOUR = SUB_PROFILE_COUNT+1;
    
    prompt_str = input(['Enter color number 1-', ...
                        num2str(NSUB_PROFILES) ' (', ...
                        num2str(DEFAULT_COLOUR), ...
                        ') : ' ],'s');
    if strcmp(prompt_str, '')
      cIndex = DEFAULT_COLOUR;  
    else
      cIndex = str2double(prompt_str);  
    end
    assert(ismember(cIndex, 1:NSUB_PROFILES),'Invalid selection for color')
    
    ls = input(['Enter a plot style (' DEFAULT_PLOTSTYLE ') : '], 's');
      
    if strcmp(ls, '')
      ls = DEFAULT_PLOTSTYLE;
    end
      
    testfig = figure;
    try
      plot([0 1], [0 1], ls);
      delete(testfig)
    catch err
      disp(err.message);
      disp('Invalid plot style')
      delete(testfig)
      throw(err)
    end
        
    CLABEL_STR = input('Enter a label ('''') : ', 's');
    if strcmp(CLABEL_STR, '')
      CLABEL_STR = '';
    end
    
    CLABEL_HA = input( ...
                'Horizontal alignment of label (left)/center/right : ', ...
                's');
      
    if strcmp(CLABEL_HA,'')
      CLABEL_HA = 'left';
    end
              
    testfig = figure;

    try
      text(xd(i), yd(i), 'hello world', ...
           'horizontalalignment', CLABEL_HA);
      delete(testfig)
    catch err
      disp(err.message);
      disp('Invalid horizontal alignment');
      delete(testfig)
      throw(err)
    end

    CLABEL_VA = input( ...
     'Vertical alignment of label top/cap/(middle)/baseline/bottom : ', ...
     's');

    if strcmp(CLABEL_VA,'')
      CLABEL_VA = 'middle';
    end
   
    testfig = figure;

    try
      text(xd(i), yd(i), 'hello world', ...
           'verticalalignment', CLABEL_VA);
      delete(testfig)
    catch err
      disp(err.message);
      disp('Invalid vertical alignment');
      delete(testfig)
      throw(err)
    end
    
    figure(src)
    
    % probably need to enforce axis about now!
    
    hold on
    plot(xd(i), yd(i), '.', 'color', [0 0 0]); %, ...
       % 'markersize', 1, 'linewidth', 2); % handl = ...
    axis(AXES_LIMITS);
    text(xd(i), yd(i), ['' CLABEL_STR ''], ...
         'verticalalignment', CLABEL_VA, ...
         'horizontalalignment', CLABEL_HA)
    drawnow
    hold off
    
    figure(branch_h_fh)
    
    hold on
    plot(xd(i), yd(i), '.', 'color', [0 0 0]); %, ...
      %  'markersize', 1, 'linewidth', 2); % handl = ...
    axis(AXES_LIMITS);
    text(xd(i), yd(i), ['$' CLABEL_STR '$'], ...
         'verticalalignment', CLABEL_VA, ...
         'horizontalalignment', CLABEL_HA)
    drawnow
    hold off
    
    figure(branch_E_fh)
    
    hold on
    plot(xd_E(i), yd_E(i), '.', 'color', [0 0 0]); %, ...
      %  'markersize', 1, 'linewidth', 2); % handl = ...
    axis(AXES_LIMITS_E);
    text(xd_E(i), yd_E(i), ['$' CLABEL_STR '$'], ...
         'verticalalignment', CLABEL_VA, ...
         'horizontalalignment', CLABEL_HA)
    drawnow
    hold off
        
    SUB_PROFILE_COUNT = SUB_PROFILE_COUNT + 1;
    ALL_PROFILE_COUNT = ALL_PROFILE_COUNT + 1;

    else
      error('Can''t find a near enough profile')
    end % waste < ...

    errmsg1=['Not gaussian pressure data, ' ...
             'use different callback function.'];

    assert(isfield(in_d, 'F'), errmsg1);
    assert(isfield(in_d, 'A'), errmsg1);
    assert(isfield(in_d, 'B'), errmsg1);
    
    % generate the name for the streamline file
   
    assert(isfield(in_d, 'phi_s_mid'), 'Missing data in m/mat-file');
    assert(isfield(in_d, 'x_s_mid'),   'Missing data in m/mat-file');
    assert(isfield(in_d, 'u_s_mid'),   'Missing data in m/mat-file');
    assert(isfield(in_d, 'v_s_mid'),   'Missing data in m/mat-file');
    
    assert(isfield(in_d, 'y_s_mid') || isfield(in_d, 'eta_s_mid'), ...
           'Missing data in m/mat-file');


    prof_phi = in_d.phi_s_mid{i};
    prof_x = in_d.x_s_mid{i};

    if isfield(in_d, 'y_s_mid')
      prof_y = in_d.y_s_mid{i};
    else
      prof_y = 1 + in_d.eta_s_mid{i};
    end

    vel_u = in_d.u_s_mid{i};
    vel_v = in_d.v_s_mid{i};
    
    prof_phi_sub = in_d.phi_sub{i};
    
    if prof_phi(1) >= 0
      % Solutions are half-waves, need to create full wave profile.
      N = length(prof_x);
      j = [N:-1:1,1:N];

      prof_phi_sub = [-in_d.phi_sub{i}(end:-1:2), in_d.phi_sub{i}];
      prof_phi = prof_phi(j);
      prof_phi(1:N) = -prof_phi(1:N);
      prof_x = prof_x(j);
      prof_x(1:N) = -prof_x(1:N);
      prof_y = prof_y(j);
      vel_u = vel_u(j);
      vel_v = vel_v(j);
      vel_v(1:N) = -vel_v(1:N);
    end
    
    disp(' ');
    disp('--------------------------------- Summary of chosen parameters');
    disp(' ');
    disp(['Froude : ' num2str(in_d.F{i})]);
    disp(['A      : ' num2str(in_d.A{i})]);
    disp(['B      : ' num2str(in_d.B{i})]);
    
    disp(' ');
    disp(['Colour number   : ' num2str(cIndex)]);
    disp(['Plot style      : ' ls]);
    

    prof_theta = atan2(vel_v, vel_u);
    %prof_tau = 0.5*log(vel_u.^2 + vel_v.^2);
    
    c = blue_black(NSUB_PROFILES);

    figure(prof_fh)
    
    hold on
    plot(prof_x,prof_y, ls, 'color', c(cIndex,:));
    
    % should add an option about printing streamlines.
    %for j = 1:length(psi_vec)
    %plot(int_X{j},psi_vec(j) + int_Y{j}, '--', 'color', [0 0 0]);
    %end
    hold off
    
    if PLOT_TOPOGRAPHY    
    % Plot the topography
      hold on
      if exist('int_xb','var')
        plot(int_xb, int_yb, ls, 'color', c(cIndex,:));
      end
      hold off
    end
    
    figure(theta_fh)
    
    hold on
    plot(prof_x,prof_theta, ls, 'color', c(cIndex,:));
    hold off
    
    if SUB_PROFILE_COUNT == NSUB_PROFILES
     % Print the output.
      disp('Sub overlay complete.')
      
      figure(prof_fh)

      % Need to add code to change ticklabels to correct font.
      axis tight
      
      ylabel('y', 'rotation', 0);
      xlabel('x');

      grid on
      box on
      
      stretch(prof_fh, [1 PROFPLOT_STRETCH]);    

      if PRINT_GRAPHS_TO_FILE
      cd('../output/')
      
      matlab2tikz(PROF_PLOTFILE_NAME, ...
                  'width','\figurewidth', ...
                  'parseStrings', false, ...
                  'figurehandle', prof_fh);
      cd('../matlab/')
      end
      
      figure(theta_fh)

      % Need to add code to change ticklabels to correct font.
      axis tight
      
      ylabel('th', 'rotation', 0);
      xlabel('x');

      grid on
      box on
      
      if PRINT_GRAPHS_TO_FILE
      cd('../output/')
      
      matlab2tikz(THETA_PLOTFILE_NAME, ...
                  'width','\figurewidth', ...
                  'parseStrings', false, ...
                  'figurehandle', theta_fh);
      cd('../matlab/')
      end
      
      SOLN_BRANCH = [];
      PLOT_TOPOGRAPHY = [];
      PROF_PLOTFILE_NAME = [];
      THETA_PLOTFILE_NAME = [];
      CB_M2TIKZ_STATE = 'CHOOSE_BRANCH';
    end
    
    if ALL_PROFILE_COUNT == total_profiles
      
      figure(branch_h_fh)

      axis(AXES_LIMITS)

      ylabel('n', 'rotation', 0);
      xlabel('w');

      grid on
      box on

      if PRINT_GRAPHS_TO_FILE
        cd('../output/')
        file_name = sprintf('%s_pressure_A%s_height_summary', ...
                            DATE_STR, astr);
        matlab2tikz(file_name, 'width','\figurewidth', ...
                               'parseStrings', false, ...
                               'figurehandle', branch_h_fh);
        cd('../matlab/')

      end
      
      figure(branch_E_fh)

      axis(AXES_LIMITS_E)

      ylabel('E', 'rotation', 0);
      xlabel('w');

      grid on
      box on

      if PRINT_GRAPHS_TO_FILE
        cd('../output/')
        file_name = sprintf('%s_pressure_A%s_energy_summary', ...
                            DATE_STR, astr);
        matlab2tikz(file_name, 'width','\figurewidth', ...
                               'parseStrings', false, ...
                               'figurehandle', branch_E_fh);
        cd('../matlab/')

      end
      
      SOLN_BRANCH = [];
      PLOT_TOPOGRAPHY = [];
      PROF_PLOTFILE_NAME = [];
      THETA_PLOTFILE_NAME = [];
      CB_M2TIKZ_STATE = 'DONE';
    end
   
  case 'DONE'
    % do nothing? add option to restart! would need to record a few things
    % tho :( related to branch plots (e.g. dots and labels etc)
    disp('Callback/summary printouts done')
  end
  
end