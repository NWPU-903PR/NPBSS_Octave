% ���ݳ���list �� �ο����� ��ȡÿ������ �� �������
% ���룺
% �����
% ���ӣ�
function [rightSeq, simulatedSeq, qv] = simulatePacBio_CCS(lengthList, genome, length_position, qv_ccs, insertMean, deleteMean, subsituteMean)
%qvNumPorpation = [0.000773775156436566,0.000902917046245173,0.00156315495789167,0.00686389144332743,0.0135195415893385,0.0144595869289036,0.0137293971602775,0.0109786749073541,0.0109915890963350,0.0116453699134911,0.0123432742096651,0.0122098275901962,0.0120661572377841,0.0125057777543409,0.0133371286699838,0.0160437274438892,0.0193820452954417,0.0250922691898122,0.0196010484169088,0.0159032856387223,0.0183101676100302,0.0218368173842203,0.0252547727344880,0.0466616552438872,0.0156955824326135,0.0161045317503407,0.0166625323325554,0.0205889838739446,0.0257030027103654,0.0333369026716767,0.0226321161889583,0.0214746820015486,0.0208989244094853,0.0186717649014943,0.0205184939257574,0.0544962632256093,0.0141071371879676,0.0106660439157758,0.00977926960575672,0.0109549988942226,0.0130826115288194,0.0159409520232498,0.0167287175510823,0.0207574064219033,0.0173711984528802,0.0157074204391793,0.0159178141013258,0.0507516865123672,0.00900926108777290,0.00731158332799726,0.00502846233442261,0.00568708597244650,0.00682138223793209,0.00748108205837106,0.00576887583599195,0.00606913072979696,0.00513392821109964,0.00427997746474023,0.00458077044975277,0.0114973948314187,0.00333885594276001,0.00313653364872653,0.00250858120953218,0.00252741440179593,0.00296111591506984,0.00308649116642569,0.00347768347430426,0.00363480610690473,0.00300846794133299,0.00248060046674031,0.00252472394575825,0.00697958105294764,0.00182574346716917,0.00186556221652683,0.00126128579046406,0.00135007083970747,0.00166108755766320,0.00120855285212554,0.00120532430488033,0.00105196831073261,0.000920674056093856,0.000831350915642903,0.000853950746359409,0.00173964887396343,0.000811441540964077,0.000540243572366003,0.000513877103196746,0.000538629298743396,0.000677456830287647,0.000628490530401884,0.000731804042248769,0.000715123214815158,0.000458991800028088,0.00575004264372820];
%qv_is = 33:126;
%errorNumProporationEachQv = [0.000501187233627273,0.000398107170553497,0.000316227766016838,0.000251188643150958,0.000199526231496888,0.000158489319246111,0.000125892541179417,0.000100000000000000,7.94328234724282e-05,6.30957344480193e-05,5.01187233627272e-05,3.98107170553497e-05,3.16227766016838e-05,2.51188643150958e-05,1.99526231496888e-05,1.58489319246111e-05,1.25892541179417e-05,1.00000000000000e-05,7.94328234724282e-06,6.30957344480193e-06,5.01187233627273e-06,3.98107170553497e-06,3.16227766016838e-06,2.51188643150958e-06,1.99526231496888e-06,1.58489319246111e-06,1.25892541179417e-06,1.00000000000000e-06,7.94328234724282e-07,6.30957344480193e-07,5.01187233627273e-07,3.98107170553497e-07,3.16227766016838e-07,2.51188643150958e-07,1.99526231496888e-07,1.58489319246111e-07,1.25892541179417e-07,1.00000000000000e-07,7.94328234724282e-08,6.30957344480193e-08,5.01187233627273e-08,3.98107170553497e-08,3.16227766016838e-08,2.51188643150958e-08,1.99526231496888e-08,1.58489319246111e-08,1.25892541179417e-08,1.00000000000000e-08,7.94328234724282e-09,6.30957344480194e-09,5.01187233627272e-09,3.98107170553497e-09,3.16227766016838e-09,2.51188643150958e-09,1.99526231496888e-09,1.58489319246111e-09,1.25892541179417e-09,1.00000000000000e-09,7.94328234724282e-10,6.30957344480194e-10,5.01187233627271e-10,3.98107170553497e-10,3.16227766016838e-10,2.51188643150958e-10,1.99526231496888e-10,1.58489319246111e-10,1.25892541179417e-10,1.00000000000000e-10,7.94328234724282e-11,6.30957344480194e-11,5.01187233627272e-11,3.98107170553497e-11,3.16227766016838e-11,2.51188643150958e-11,1.99526231496888e-11,1.58489319246111e-11,1.25892541179417e-11,1.00000000000000e-11,7.94328234724282e-12,6.30957344480194e-12,5.01187233627272e-12,3.98107170553497e-12,3.16227766016838e-12,2.51188643150958e-12,1.99526231496888e-12,1.58489319246111e-12,1.25892541179417e-12,1.00000000000000e-12,7.94328234724282e-13,6.30957344480194e-13,5.01187233627272e-13,3.98107170553497e-13,3.16227766016838e-13,2.51188643150958e-13];
gen = upper(genome); % ���ַ�ȫ��ת��Ϊ��д
genomeLength = length(gen);
num = length(lengthList);
rightSeq = cell(num,1);
simulatedSeq = cell(num,1);
qv = cell(num,1);

insertStd  = 0.0165;

deleteStd  = 0.01084;

subsituteStd = 0.01585;


%totalErrors  = normrnd(errorMean, errorStd, num);
insertErrors = abs(normrnd(insertMean, insertStd, num, 1));
deleteErroes = abs(normrnd(deleteMean, deleteStd, num, 1));
subsituteErrors = abs(normrnd(subsituteMean, subsituteStd, num, 1));
for i = 1 : num
    rr = unifrnd(1,genomeLength - lengthList(i));
    rr = ceil(rr);
    ss = gen(rr : rr + lengthList(i) - 1);
    s = nt2int(ss) ;
    rightSeq{i} = ss;
    insertNum = ceil(insertErrors(i) * lengthList(i));
    deleteNum = ceil(deleteErroes(i) * lengthList(i));
    subsituteNum = ceil(subsituteErrors(i) * lengthList(i));
    errorNum = insertNum + deleteNum + subsituteNum;
%     
%     qvFirst = 94 * ones(1, lengthList(i));
%     qvNum = lengthList(i) * qvNumPorpation;
%     qvErrorNum = qvNum .* errorNumProporationEachQv;
%     qvNumPorpationActually = qvErrorNum/sum(qvErrorNum);
    position = randperm(lengthList(i), errorNum);
    insertPosition = position(1:insertNum);
    deletePosition = position(insertNum + 1 : insertNum + deleteNum);
    subsitutePosition = position(insertNum + deleteNum + 1 : end);
    %subsituteRand = randint(1,subsituteNum,[1, 4]);
    % for subsitution
    s(subsitutePosition) = s(subsitutePosition) + 1;
    s(find(s>4)) = s(find(s>4)) - 4;
    
    % for delete
    s(deletePosition) = 0;
    
    % for insert
    insertLetters = randi([1, 4],1,insertNum);
    insertPosition = sort(insertPosition);
    
    qvis = qv_ccs{length_position(i)};
    %assignQVs(errorNum, position, qvFirst, qvNumPorpationActually);
    %afterQV = qvis(1 : insertPosition(1));
    afters = s(1 : insertPosition(1));
    afters = [afters, insertLetters(1)];
    %afterQV = [afterQV, 0];
    index = insertPosition(1) + 1;
    for j = 2 : insertNum
        afters = [afters, s(index : insertPosition(j))];
        afters = [afters, insertLetters(j)];
        %fprintf('j = %d, len(afters) - length(sub) = %d\n',j,length(afters) - insertPosition(j));
        %afterQV = [afterQV, qvis(index : insertPosition(j))];
        %afterQV = [afterQV, 0];
        
        index = insertPosition(j) + 1;
    end
    afters = [afters, s(index : end)];
    %afterQV = [afterQV, qvis(index : end)];
    %afterQV(deletePosition) = -1;
    %finalQV = afterQV(find(afterQV>=0));
    
    finals = afters(find(afters>0));
    qvff = [qvis,qvis];
    finalQV = qvff(1:length(finals));
    simulate = int2nt(finals,'Case','upper');
    simulatedSeq{i} = simulate;
    qv{i} = finalQV;
end
end