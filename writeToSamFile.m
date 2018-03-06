% write to fasta file
% input:
% output:
% example:
function writeToSamFile(samFormat,qv, fileName, genomeName, genomeLength)
outfilename = fileName;
fid=fopen(outfilename,'w');
firstt = ['@SQ\tSN:', genomeName, '\t', 'LN:', num2str(genomeLength), '\n'];
fprintf(fid,firstt);
for i=1:length(samFormat)
      fprintf(fid,samFormat{i});
      q = qv{i};
      phred = q + 33;
      strr = char(phred);
      fprintf(fid,'%s\t\n',strr);
% fprintf(fid,'%s\n',Clones(i,:));
end
fclose(fid);
end