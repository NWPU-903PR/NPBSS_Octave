function x = randsample(mat,n)
  if n <= length(mat)
    x=zeros(n,1);
    indexx = randperm(length(mat));
    x = mat(indexx(1:n));
  else
    for i = 1:n
      indexx = randi(length(mat),1);
      x(i) = mat(indexx);
    end
  end
end