function [a] = unforced_integral_check_b(in_d, i, exact_d)

  f1 = @(x, fs, fs_exact)(abs(ppval(fs,x) - ppval(fs_exact,x)));
  
  c_vals = cell2mat(exact_d.cluster_val);
  
  % construct data for delta integrand
  delta_int = zeros(size(c_vals));
  for k = 1:length(c_vals)
    delta_int(k) = quad(@(x)f1(x, ...
                               in_d.eta_s_pp{i}, ...
                               exact_d.eta_s_pp{k}), ...
                        0, exact_d.x_s_mid{k}(end));
  end
  
  % construct spline for delta integrand
  delta_int_pp = spline(c_vals, delta_int);
  
  f2 = @(x, d_int)(ppval(d_int, x));
  a = quad(@(x)f2(x, delta_int_pp), ...
           min(c_vals), max(c_vals)) / (max(c_vals) - min(c_vals));

end