% Konstantinos Letros 8851
% Multimedia Systems Project
% Inverse Run Length Encoding

function qBlock = irunLength(runSymbols, DCpred)

% Part 1
vec = [];

% Inverse RLE
for i=1:length(runSymbols)-1
    vec = [vec;zeros(runSymbols(i,1),1)];
    vec(end+1,1)=runSymbols(i,2);
end

% Adding EOB zeros
vec = [vec;zeros(64-length(vec),1)];

% Part 2
global zigZagTable
Tables;

qBlock = zeros(sqrt(length(vec)));

% Place values back, in zig-zag form
for i = 1:length(vec)
    [r,c]=find(zigZagTable==i);
    qBlock(r,c) = vec(i);
end

qBlock(1,1) = qBlock(1,1) + DCpred;

end

%% Huffman Encoding
function huffStream = huffEnc(runSymbols)

global Nx
binaryStream = [];

for k = 1:size(runSymbols,1)
    
    % Find the Category of the non-zero DCT coeffs
    category = 0;
    
    if runSymbols(k,2)~=0
        while category<11
            if 2^(category-1)<=abs(runSymbols(k,2)) && abs(runSymbols(k,2))<=(2^category)-1
                break
            end
            category = category + 1;
        end
    else
        category = 0;
    end
    
    % Convert non-zero DCT Coeffs to binary
    
    if runSymbols(k,2)<0 % 1's Complement
        magnitude = double(~de2bi(-runSymbols(k,2),'left-msb'));
    else
        magnitude = de2bi(runSymbols(k,2),'left-msb');
    end
    
    
    if k==1
        % Find Index of Category in DC table
        huffmanInd = category + 1;
        DCvector=[Nx.DCTable{huffmanInd}-'0',magnitude];
        binaryStream = [binaryStream,DCvector];
    else
        % Find Index of Run/Category in AC table
        huffmanInd = 10*runSymbols(k,1)+(category+1);
        if runSymbols(k,1)==15 % index 152 --> ZRL
            huffmanInd = huffmanInd + 1;
        elseif isequal(runSymbols(k,:),[0 0]) % EOB
            magnitude = [];
        end
        ACvector = [Nx.ACTable{huffmanInd}-'0',magnitude];
        binaryStream = [binaryStream,ACvector];
    end
    
end

% Zero padding to create bytes
while (mod(length(binaryStream),8)) ~= 0
    binaryStream(end+1)=0;
end

% Convert Stream to bytes
i = 0;
k = 1;
huffStream = zeros(1,length(binaryStream)/8);
while i < length(binaryStream)
    huffStream(k) = uint8(bi2de(binaryStream(i+1:i+8),'left-msb'));
    i = i + 8;
    k = k + 1;
end

end

%% Huffman Decoding
function runSymbols = huffDec(huffStream)

global DCTable
global ACTable

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
while length(find(DCTable==strStream))~=1
    stop = stop + 1;
    strStream = sprintf('%d',binaryStream(bit:stop));
end

category = find(DCTable==strStream)-1;
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
    while length(find(ACTable==strStream))~=1
        stop = stop + 1;
        strStream = sprintf('%d',binaryStream(bit:stop));
    end
    
    if isequal(strStream , ACTable{1}) % EOB
        runSymbols(k,:)= [0,0];
        break;
    end
    
    indexAC = find(ACTable==strStream)-1;
    precZeros = fix(indexAC/10);
    category = mod(indexAC,10);
    if indexAC > 151
        category = category - 1;
    end
    
    bit = stop + 1;
    stop = stop + category;
    
    if isequal(strStream , ACTable{152}) % ZRL
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