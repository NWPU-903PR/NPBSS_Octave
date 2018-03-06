function x = nt2int(ss)
  len = length(ss);
  x = ones(1, len);
  for i = 1:len
    if ss(i) == 'C'
      x(i) = 2;
    elseif ss(i) == 'G'
      x(i) = 3;
    elseif ss(i) == 'T'
      x(i) = 4;
    end
  end
end  