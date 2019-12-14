clc
clear

zzBlockVec = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% zzBlockVec = zeros(64,1);

zzBlockVec = [zzBlockVec ; NaN]; % Adding delimeters


lastNonZeroInd = 1;
runSymbols(1,1) = 0;
runSymbols(1,2) = zzBlockVec(1);
k = 2;
i = 2;

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


sum(vec == zzBlockVec(1:end-1))/length(vec)