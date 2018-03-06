function x = int2nt(ss)
  len = length(ss);
  x = '';
  for i = 1:len
    if ss(i) == 1
      x = [x, 'A'];
    elseif ss(i) == 2
      x = [x, 'C'];
    elseif ss(i) == 3
      x = [x, 'G'];
    elseif ss(i) == 4
      x = [x, 'T'];
    end
  end
end  