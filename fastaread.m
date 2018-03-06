function [head, genome] = fastaread(referenceFile)
  fidin=fopen(referenceFile); 
  head = fgetl(fidin);
  head = head(2:end);
  genome = fgetl(fidin);
  fclose(fidin);
  end