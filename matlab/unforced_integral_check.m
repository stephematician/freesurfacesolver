function [a] = unforced_integral_check(in_d, i)

  F = in_d.F{i};
  A = in_d.A{i};
  B = in_d.B{i};

  eta = in_d.eta_s_mid{i} + 0.5*(F^2)*(1 - (in_d.u_s_mid{i}(end)^2 + in_d.v_s_mid{i}(end)^2));
  
  free_surf_pp = spline(in_d.x_s_mid{i},eta);
  dfree_surf_pp = fnder(free_surf_pp, 1);
  
  f1 = @(fs, x, F)((3/(2*F^2))*ppval(fs, x).^2 + (F^(-2) - 1)*ppval(fs, x));
  
  function [b] = f2(fs, dfs, X, F, A, B, x_end)
    psub = @(x, A, B)(A * B * exp(-(B*x).^2) / sqrt(pi));
    f2sub = @(dfs, x, F, A, B)(ppval(dfs, x) .* psub(x, A, B));

    b_0 = zeros(size(X));
    j = 1;
    for x = X
      b_0(j) = (2*(1+ppval(fs, x)).*psub(x, A, B) - ...
                 quad(@(y)f2sub(dfs, y, F, A, B), x, x_end)) / (F^2);
      j = j + 1;
    end
    b = b_0;
  end

  a1 = 2*quad(@(x)f1(free_surf_pp, x, F), 0, in_d.x_s_mid{i}(end));
  a2 = quad(@(x)f2(free_surf_pp, dfree_surf_pp, x, F, A, B, in_d.x_s_mid{i}(end)), 0, in_d.x_s_mid{i}(end));

  a = a1 + a2;

end