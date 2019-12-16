% Konstantinos Letros 8851
% Multimedia Systems Project
% Huffman Decoding

function runSymbols = huffDec(huffStream,type)

global DCCategoryCode
global ACCategoryCode

% Convert bytes to binary stream
binaryStream = [];
for i = 1:length(huffStream)
    binaryStream = [binaryStream,de2bi(huffStream(i),8,'left-msb')];
end

% Initialize
k=1;
bit = 1;
stop = 2;

% Convert binary vector to string
strStream = sprintf('%d',binaryStream(bit:stop));

% DC Difference Coefficient
% Repeat until unique reference (DC Table)
while length(find(DCCategoryCode{type}==strStream))~=1
    stop = stop + 1;
    strStream = sprintf('%d',binaryStream(bit:stop));
end

category = find(DCCategoryCode{type}==strStream)-1;
bit = stop + 1;
stop = stop + category;

if category == 0 % DC = 0
    DCdiff = 0;
    stop = stop + 1;
else
    % 0 -> Negative Sign
    if binaryStream(bit)==0
        DCdiff = -bi2de(~binaryStream(bit:stop),'left-msb');
    else
        DCdiff = bi2de(binaryStream(bit:stop),'left-msb');
    end
end
% Add DC to runSymbols
runSymbols(k,:)= [0 , DCdiff];

stop = stop + 1;
bit = stop;
k = k + 1;

% AC Coefficients
while bit < length(binaryStream)
    
    % Convert binary vector to string
    strStream = sprintf('%d',binaryStream(bit:stop));
    
    % Repeat until unique reference (ACTable)
    while length(find(ACCategoryCode{type}==strStream))~=1
        stop = stop + 1;
        strStream = sprintf('%d',binaryStream(bit:stop));
    end
    
    if isequal(strStream , ACCategoryCode{type}{1}) % EOB
        runSymbols(k,:)= [0,0];
        break;
    end
    
    indexAC = find(ACCategoryCode{type}==strStream)-1;
    precZeros = fix(indexAC/10);
    category = mod(indexAC,10);
    if indexAC > 151
        category = category - 1;
    end
    
    bit = stop + 1;
    stop = stop + category;
    
    if isequal(strStream , ACCategoryCode{type}{152}) % ZRL
        runSymbols(k,:)= [15,0];
        k = k + 1;
        stop = stop + 1;
        bit = stop;
        continue;
    end
    
    % 0 -> Negative Sign
    if binaryStream(bit)==0
        AC = -bi2de(~binaryStream(bit:stop),'left-msb');
    else
        AC = bi2de(binaryStream(bit:stop),'left-msb');
    end
    
    runSymbols(k,:)= [precZeros,AC];
    k = k + 1;
    stop = stop + 1;
    bit = stop;
end

end