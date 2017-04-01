function [ a ] = stretch( f, factor )
  m = gcf();

  figure(f)

  c_axis = axis();
  
  if length(c_axis) == 6
    w = diff(c_axis);
    w = w(1:2:5);

    c_axis(1:2:5) = c_axis(1:2:5) - (w .* (factor - 1) / 2);
    c_axis(2:2:6) = c_axis(2:2:6) + (w .* (factor - 1) / 2);

  elseif length(c_axis) == 4
    w = diff(c_axis);
    w = w([1 3]);

    c_axis([1 3]) = c_axis([1 3]) - (w .* (factor - 1) / 2);
    c_axis([2 4]) = c_axis([2 4]) + (w .* (factor - 1) / 2);  
  end

  axis(c_axis);

  figure(m);
end
