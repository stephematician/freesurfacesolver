function [x,y,xval,yval,h] = batch_plot_m2tikz(input, f1, f2, varargin)
% Plots data extracted from free-surface using interpolating or smoothing
% splines.
%
% Stephen Wade
% 12/12/2013

% To do
% Need to fix up the number of points this thing plots, it think it gets a
% bit yuk, especially when the data is so large now.

% History
% 12/12/2013 Adapted from batch_plot_sp.m

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

  if isfield(input, 'x') 
    n_surface = length(input.x);
  elseif isfield(input, 'x_mid')
    n_surface = length(input.x_mid);
  elseif isfield(input, 'x_s_mid')
    n_surface = length(input.x_s_mid);
  end
  
  xval = zeros(1, n_surface);
  yval = zeros(1, n_surface);
  
  for i = 1 : n_surface
    xval(i) = f1(input, i);
    yval(i) = f2(input, i);
  end
  
  ref_ds = sqrt((max(xval) - min(xval)).^2 + (max(yval) - min(yval)).^2);
  
  % I wish matlab did regexp's here.
  xc = {};
  yc = {};
  i = 1;
  k = find(~isnan(xval),1,'first');
  while ~isempty(k)
    j = find(isnan(xval(k+1:end)),1,'first');
    if ~isempty(j)
      xc{i} = xval(k:k+(j-1));
      yc{i} = yval(k:k+(j-1));
    else
      xc{i} = xval(k:end);
      yc{i} = yval(k:end);
      break
    end
    k = (k + j) + find(~isnan(xval(k+j+1:end)),1,'first');
    i = i + 1;
  end
  
  h = zeros(size(xc));
  xs = cell(size(xc));
  ys = cell(size(xc));
  
 
  for i = 1:length(xc)
    % Need a kind of arclength variable
    ds = [0, cumsum(sqrt(diff(xc{i}).^2 + diff(yc{i}).^2))];
    if ~isOctave
      fitf = fit(ds(:), xc{i}(:), 'smoothingspline','SmoothingParam',1-4e-11);
      xs{i} = feval(fitf,linspace(min(ds),max(ds), 200*max(ds)/ref_ds));
      fitf = fit(ds(:), yc{i}(:), 'smoothingspline','SmoothingParam',1-4e-11);
      ys{i} = feval(fitf,linspace(min(ds),max(ds), 200*max(ds)/ref_ds));
    else
      xs{i} = spline(ds,xc{i},linspace(min(ds),max(ds),200*max(ds)/ref_ds));
      ys{i} = spline(ds,yc{i},linspace(min(ds),max(ds),200*max(ds)/ref_ds));
    end
    %xs{i} = csaps(ds,xc{i},[],linspace(min(ds),max(ds),length(ds)*20)); alternative octave option
    %ys{i} = csaps(ds,yc{i},[],linspace(min(ds),max(ds),length(ds)*20));

    h(i) = plot(xs{i}, ys{i}, varargin{:});
  end

  x = cell(size(xs));
  y = cell(size(xs));
  
  for i = 1:length(xs)
    x{i} = xs{i};
    y{i} = ys{i};
  end

end
