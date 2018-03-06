function [head, seq, qv] = fastqread(referenceFile)
  fidin=fopen(referenceFile); 
  head = {};
  seq = {};
  qv = {};
  n = 1;
    while  ~feof(fidin)
      head{n} = fgetl(fidin);
      seq{n}  = fgetl(fidin);
      qv{n}   = fgetl(fidin);
      qv{n}   = fgetl(fidin);
      n = n + 1;
  
    end
    fclose(fidin);
end