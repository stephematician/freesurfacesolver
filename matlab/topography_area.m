function [a] = topography_area(in_d, i)

  N = length(in_d.x_s_mid{i});
  j = [N:-1:1,1:N];

  prof_phi = in_d.phi_s_mid{i};


  prof_phi_sub = [-in_d.phi_sub{i}(end:-1:2), in_d.phi_sub{i}];
  prof_phi = prof_phi(j);
  prof_phi(1:N) = -prof_phi(1:N);

  vel_u = exp(in_d.tau_s_mid{i}) .* cos(in_d.theta_s_mid{i});
  vel_v = exp(in_d.tau_s_mid{i}) .* sin(in_d.theta_s_mid{i});

  vel_u = vel_u(j);
  vel_v = vel_v(j);
  vel_v(1:N) = -vel_v(1:N);

  prof_theta = atan2(vel_v, vel_u);

  [ int_xb, ...
    int_yb, ...
    ~, ...
    ~, ...
    ~] = compute_htopography_floor(prof_phi, ...
                               prof_phi_sub, ...
                                  in_d.A{i}, ...
                                  in_d.B{i}, ...
                                  in_d.L{i}, ...
                              in_d.phi_c{i}, ...
                                 prof_theta);

  a = trapz(int_xb, int_yb);

end