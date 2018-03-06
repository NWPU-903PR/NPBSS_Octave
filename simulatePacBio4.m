% QV table and model can update
function [rightSeq, simulatedSeq, qv, samFormat] = simulatePacBio4(lengthList, genome, genomeName,ins, del, sub, samFlag, qvNumPorpation, errorNumProporationEachQv)
% qvNumPorpation = [0.000276726086402811,0.00832439676688714,0.0253280741491020,0.0323037535959514,0.0306580764541609,...
%                   0.0359438422593862,0.0405823269033696,0.0479061456202096,0.0606710538165784,0.0689325673079174,...
%                   0.0724403432117603,0.0895914094701023,0.121460336125986,0.206006999285471,0.159571766877217];
%errorNumProporationEachQv = [0.3634,0.3421,0.2868,0.1766,0.1482,0.1127,0.0978,0.0869,0.0686,0.0584,0.0555,0.0486,0.0426,0.0307,0.0207];
qvTable = qvNumPorpation;
model = errorNumProporationEachQv;
gen = upper(genome); % 将字符全部转化为大写
qv_uniform = zeros(100000,1);
randomp = randperm(100000);
aaa = floor(100000*qvNumPorpation);
haha = 1;
for ii = 1 : length(aaa)
    indexx = randomp(haha : haha + aaa(ii));
    qv_uniform(indexx) = ii - 1;
    haha = haha + aaa(ii);
end
genomeLength = length(gen);
num = length(lengthList);
rightSeq = cell(num,1);
sampleIndex = randsample([1 2 3 4],num,1);
yyy=@(x)[45-60000/(x+7500)];%[44-28000/(x+3500)];
simulatedSeq = cell(num,1);
samFormat = cell(num,1);
qv = cell(num,1);
insertMean = ins;
insertStd  = 0.0212;
deleteMean = del;  deleteStd  = 0.02761;
subsituteMean = sub; subsituteStd = 0.01857;
%totalErrors  = normrnd(errorMean, errorStd, num);
insertErrors = abs(normrnd(insertMean, insertStd, num, 1));
deleteErroes = abs(normrnd(deleteMean, deleteStd, num, 1));
subsituteErrors = abs(normrnd(subsituteMean, subsituteStd, num, 1));
%MyPar = parpool(2); 
for i = 1 : num
    rr = unifrnd(1,genomeLength - lengthList(i));
    rr = ceil(rr);
    ss = gen(rr : rr + lengthList(i) - 1);
    s = nt2int(ss) ;
    rightSeq{i} = ss;
    insertNum = ceil(insertErrors(i) * lengthList(i));
    deleteNum = ceil(deleteErroes(i) * lengthList(i));
    subsituteNum = ceil(subsituteErrors(i) * lengthList(i));
    insertPosition = randperm(lengthList(i) - 1, insertNum) + 1;
    insertLetters = randi([1, 4],1,insertNum);
    insertPosition = sort(insertPosition);
    insertPosition2 = insertPosition;%存放插入后在序列中的insert位置
    afters = s(1 : insertPosition(1));
    afters = [afters, insertLetters(1)];
    insertPosition2(1) = insertPosition2(1) + 1;
    index = insertPosition(1) + 1;
    for j = 2 : insertNum
        afters = [afters, s(index : insertPosition(j))];
        afters = [afters, insertLetters(j)];
        %fprintf('j = %d, len(afters) - length(sub) = %d\n',j,length(afters) - insertPosition(j));
        index = insertPosition(j) + 1;
        insertPosition2(j) = insertPosition(j) + j;
    end
    afters = [afters, s(index : end)];
    allPositionNum = lengthList(i) + insertNum;
    positionRemain = setdiff(2:allPositionNum - 1,insertPosition2);
    position = randsample(positionRemain, subsituteNum + deleteNum);
    deletePosition = position(1: deleteNum);
    subsitutePosition = position(deleteNum + 1 : end);
    matchPosition = setdiff(1:allPositionNum,insertPosition2);
    matchPosition = setdiff(matchPosition,deletePosition);
    matchPosition = setdiff(matchPosition,subsitutePosition);
    errorIndex = zeros(1, allPositionNum);
    errorIndex(matchPosition)    = 21; % match的位置显示21,,，方便Sam格式
    errorIndex(insertPosition2)  = 22; % insert位置显示22
    errorIndex(deletePosition)   = 23; % delete位置显示23
    errorIndex(subsitutePosition)= 24; % subsitute位置显示24
    % for subsitution
    afters(subsitutePosition) = afters(subsitutePosition) + 1;
    afters(find(afters>4)) = afters(find(afters>4)) - 4;
    if samFlag
        samFormat_each = ['Simulated_' num2str(i) '\t' '0' '\t' genomeName '\t' num2str(rr) '\t100\t'];
        matchNum = 1;
        Inum = 0;
        Dnum = 0;
        Xnum = 0;
        CIGAR = [];
        for k = 2 : allPositionNum

            if errorIndex(k) == 21   % match
                matchNum = matchNum + 1;
                if errorIndex(k - 1) == 22
                    CIGAR=[CIGAR, num2str(Inum), 'I'];
                    Inum = 0;
                elseif errorIndex(k - 1) == 23
                    CIGAR=[CIGAR, num2str(Dnum), 'D'];
                    Dnum = 0;
                elseif errorIndex(k - 1) == 24
                    CIGAR=[CIGAR, num2str(Xnum), 'X'];
                    Xnum = 0;
                end

            elseif errorIndex(k) == 22 % insert
                Inum = Inum + 1;
                if errorIndex(k - 1) == 21
                    CIGAR=[CIGAR, num2str(matchNum), '='];
                    matchNum = 0;
                elseif errorIndex(k - 1) == 23
                    CIGAR=[CIGAR, num2str(Dnum), 'D'];
                    Dnum = 0;
                elseif errorIndex(k - 1) == 24
                    CIGAR=[CIGAR, num2str(Xnum), 'X'];
                    Xnum = 0;
                end
            elseif errorIndex(k) == 23 % delete
                Dnum = Dnum + 1;
                if errorIndex(k - 1) == 21
                    CIGAR=[CIGAR, num2str(matchNum), '='];
                    matchNum = 0;
                elseif errorIndex(k - 1) == 22
                    CIGAR=[CIGAR, num2str(Inum), 'I'];
                    Inum = 0;
                elseif errorIndex(k - 1) == 24
                    CIGAR=[CIGAR, num2str(Xnum), 'X'];
                    Xnum = 0;
                end
            elseif errorIndex(k) == 24 % subsititute
                Xnum = Xnum + 1;
                if errorIndex(k - 1) == 21
                    CIGAR=[CIGAR, num2str(matchNum), '='];
                    matchNum = 0;
                elseif errorIndex(k - 1) == 22
                    CIGAR=[CIGAR, num2str(Inum), 'I'];
                    Inum = 0;
                elseif errorIndex(k - 1) == 23
                    CIGAR=[CIGAR, num2str(Dnum), 'D'];
                    Dnum = 0;                
                end
            end
        end
    
        CIGAR=[CIGAR, num2str(matchNum), '='];
        samFormat_each = [samFormat_each '\t' CIGAR '\t=\t0\t' int2str(allPositionNum - deleteNum) '\t' int2nt(afters,'Case','upper')];
    else
        samFormat_each = [];
    end
    % for delete
    afters(deletePosition) = 0;
    insertAndSubsititutePosition = union(subsitutePosition, insertPosition2);
    errorAssignedNum = 0;
    similutedNum = allPositionNum - deleteNum;
    insertAndSubError = insertNum + subsituteNum;
    % assign QV
    qvis = qv_uniform(1 : allPositionNum)';
    for m = 1 : 15
        qvErrorNum = ceil(model(m)*similutedNum*qvTable(m)/(1-model(m)));
        if errorAssignedNum + qvErrorNum > insertAndSubError
            qvErrorNum = insertAndSubError - errorAssignedNum;
        end        
        positionIndex2 = insertAndSubsititutePosition(errorAssignedNum + 1 : errorAssignedNum + qvErrorNum);
        qvis(positionIndex2) = m - 1;
        errorAssignedNum = errorAssignedNum + qvErrorNum;
        if errorAssignedNum >= insertAndSubError
            break;
        end
        
    end
    %fprintf('i = %d\n',i);
    finals  = afters(find(afters>0));
    finalQV = qvis(find(afters>0));
    simulate = int2nt(finals,'Case','upper');
    simulatedSeq{i} = simulate;
    if sampleIndex(i) == 4
        fitValue = floor(yyy(lengthList(i))) - 1;
        randomQV = randsample(fitValue:43,1) - 33;
        finalQV = qv_uniform(1:length(finalQV))';
        proportionNum = length(finalQV);
        %proportionNum = floor(length(finalQV) * randsample(1:ceil(100*((43-yyy(lengthList(i)))/7)),1)/100);
        if randomQV == 3
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/7)),1)/100);
        elseif randomQV == 4
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/6)),1)/100);
        elseif randomQV == 5
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/5)),1)/100);
        elseif randomQV == 6
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/4)),1)/100);
        elseif randomQV == 7
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/3)),1)/100);
        elseif randomQV == 8
            proportionNum = floor(length(finalQV) * randsample(0:ceil(100*((44-yyy(lengthList(i)))/2)),1)/100);
        end
        finalQV(1:proportionNum) = randomQV;
    end
    qv{i} = finalQV;
    if samFlag
        phred = finalQV + 33;
        strr = char(phred);
        samFormat_each = [samFormat_each, '\t', strr];
    end
    samFormat{i} = samFormat_each;
    
    
end
%delete(MyPar) %计算完成后关闭并行处理池
end

% assign qv to each position
function qvis = assignQVs(errorNum, errorPosition, qvFirst, eachQVProporation)
qvis = qvFirst;
eachQVNum = floor(errorNum * eachQVProporation);
n = 0;
for i = 1 : 15
   positionIs = errorPosition(n + 1 : n + eachQVNum(i));
   qvis(positionIs) = i - 1;
   n =n + eachQVNum(i);
end

end