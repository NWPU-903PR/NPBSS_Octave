% 统计CCS 序列的qv值分布，每个位置
%[heads seqs qvs]=fastqread('H:\Pacbio data\Chunlab_mock_16S\m130426_222913_42141_c100519052550000001823079909281320_s1_p0.1.subreads.fastq');

function qv_average_each_read(file)

[heads, seqs, qvs] = fastqread(file);
n_seq = length(qvs);
seq_length = zeros(1,n_seq);
average_qv = zeros(1,n_seq);
average_qv2 = zeros(1,n_seq);
qvss = qvs;
%max_length = 0;
%is_or_not = randperm(2,n_seq);
for i = 1 : n_seq
    qvss{i} = qvs{i};
    seq_length(i) = length(qvs{i});
    average_qv(i) = mean(qvss{i}+0);

end
scatter(seq_length,average_qv, 4, 'filled')
xlabel('Read Length');
ylabel('Read Average Quality');


