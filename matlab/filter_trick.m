function [a] = filter_trick(f)
  a = zeros(size(f));
 % warning('off','Octave:divide-by-zero');
  a(~f) = 0/0;
  a(f) = 1;
 % warning('on','Octave:divide-by-zero');
end
