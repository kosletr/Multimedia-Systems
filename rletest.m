clc
clear

% zzBlockVec = [121 30 0 0 3 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0]';
zzBlockVec = zeros(64,1);
length(zzBlockVec)
i = 1;
k = 1;
lastNonZeroInd = 0;

zzBlockVec = [zzBlockVec;NaN]; % Adding ending delimeter

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

% Adding EOB
runSymbols=runSymbols(1:lastNonZeroInd,:); 
runSymbols(end+1,:) = [0,0];

runSymbols


% Part 1
vec = [];

% Inverse RLE
for i=1:length(runSymbols)-1
    vec = [vec;zeros(runSymbols(i,1),1)];
    vec(end+1,1)=runSymbols(i,2);
end
vec = [vec;zeros(64-length(vec),1)];
vec