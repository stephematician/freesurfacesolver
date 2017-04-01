function [x_coarse, y_coarse] = coarse_grid_approx(x, y, tolerance, max_iterations)

  x_coarse = x(:)';
  y_coarse = y(:)';
  
  norm_xy = trapz(x,abs(y));
  
  iterations = 0;

  max_remove_ind = 2;
  
  best_error = 0;
  
  while (max_remove_ind > 1) && (iterations < max_iterations) 

  remove_inds = 2:2:(length(x_coarse)-1);
  error_vals = zeros(1, length(remove_inds));
  % okay so we look to remove every second point, and see how far we get!
  k = 1;
  for j=2:2:(length(x_coarse)-1)
    inds = [1, 2:(j-1), (j+1):(length(x_coarse)-1), length(x_coarse)];
    x_temp = x_coarse(inds);
    y_temp = y_coarse(inds);
    
    error_vals(k) = trapz(x, abs(y - interp1(x_temp, y_temp, x, 'linear')));
    k = k + 1;
  end

  % remove as much as we can in one fell swoop
  [sorted_error_vals, sorted_i] = sort(error_vals, 'ascend');
  sum_error_vals = cumsum(sorted_error_vals-best_error);
  max_remove_ind = find((best_error+sum_error_vals) >= tolerance*norm_xy, 1, 'first');
  if (isempty(max_remove_ind))
    max_remove_ind = length(sum_error_vals+1);
  end

  inds = ones(1,length(x_coarse));
  inds(remove_inds(sorted_i(1:(max_remove_ind-1)))) = 0;
  
  inds = logical(inds);
  
  x_coarse = x_coarse(inds);
  y_coarse = y_coarse(inds);
    
  best_error = trapz(x, abs(y - interp1(x_coarse, y_coarse, x, 'linear')));

  iterations = iterations + 1;
  
  end
  
end