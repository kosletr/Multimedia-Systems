% Konstantinos Letros 8851
% Multimedia Systems Project
% Inverse Run Length Encoding

function qBlock = irunLength(runSymbols, DCpred)

% Part 1
vec = [];

% Inverse RLE
for i=1:length(runSymbols)
    vec = [vec;zeros(runSymbols(i,1),1)];
    vec(end+1,1)=runSymbols(i,2);
end

% Adding EOB zeros
vec = [vec;zeros(64-length(vec),1)];

% Part 2
global zigZagTable
% Tables;

qBlock = zeros(sqrt(length(vec)));

% Place values back, in zig-zag form
for i = 1:length(vec)
    [r,c]=find(zigZagTable==i);
    qBlock(r,c) = vec(i);
end

qBlock(1,1) = qBlock(1,1) + DCpred;

end