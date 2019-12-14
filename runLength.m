% Konstantinos Letros 8851
% Multimedia Systems Project
% Run Length Encoding

function runSymbols = runLength(qBlock, DCpred)

% Part 1
global zigZagTable
Tables;

zzBlockVec = zeros(length(qBlock)^2,1);

% Get values in zig-zag form
for i = 1:length(qBlock)
    for j = 1:length(qBlock)
        zzBlockVec(zigZagTable(i,j),1) = qBlock(i,j);
    end
end
zzBlockVec(1) = zzBlockVec(1) - DCpred;

% Part 2
zzBlockVec = [zzBlockVec ; NaN]; % Adding delimeters

% DC Term
lastNonZeroInd = 1;
runSymbols(1,1) = 0;
runSymbols(1,2) = zzBlockVec(1);
k = 2;
i = 2;

% RLE
while i < length(zzBlockVec)+1
    countZeros = 0;
    while zzBlockVec(i)==0 && countZeros<15
        countZeros = countZeros + 1;
        i = i + 1;
    end
    runSymbols(k,1)=countZeros;
    runSymbols(k,2)=zzBlockVec(i);
    if zzBlockVec(i) ~= 0 && isnan(zzBlockVec(i))==false
        lastNonZeroInd = k;
    end
    k = k + 1;
    i = i + 1;
end

% Removing EOB zeros
runSymbols=runSymbols(1:lastNonZeroInd,:);
runSymbols(end+1,:) = [0,0]; % EOB

end