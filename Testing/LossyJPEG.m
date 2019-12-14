% Konstantinos Letros 8851
% Multimedia Systems Project

%% Clean the screen

clc
clear

%% Testing

originalImage = imread('flower2.jpg');
figure
imshow(originalImage)

subimg = [4,2,0];
qScale = 1;

JPEGenc = JPEGencode(originalImage,subimg,qScale);
imgREC = JPEGdecode(JPEGenc,subimg,qScale);
figure
imshow(imgREC)

% [Y,Cb,Cr] = convert2ycbcr(originalImage,subVec);
% RGBimage = convert2rgb(Y,Cb,Cr, subVec);
%
% figure
% imshow(RGBimage)

% i=randi([0,239]);
% j=randi([0,134]);
% Block = originalImage(i+1:i+8,j+1:j+8,1);
% dctBlock = blockDCT(Block);
% qBlock = quantizeJPEG(dctBlock,qTable{1},1);
% runSymbols = runlength(qBlock,0);
% huffstream = huffEnc(runSymbols,1);
%
% runSymbols = huffDec(huffstream,1);
% qblock = irunlength(runSymbols,0);
% dctblock = dequantizeJPEG(qBlock,qTable{1},1);
% block = iblockDCT(dctblock);
% block = uint8(block);
%
% sum(sum(qblock==qBlock))/64

%% Convert RGB Image to YCbCr Image
function [imageY, imageCb, imageCr] = convert2ycbcr(imageRBG, subimg)

% Transformation Matrix
Tycbcr = [.299,.587,.114;
    -.168736,-.331264,.5;
    .5,-.418688,-.081312];

newImg = zeros(size(imageRBG));
n = zeros(3,1);

% Transformation from RGB to YCbCr
for i = 1:size(imageRBG,1)
    for j = 1:size(imageRBG,2)
        n(:,1) = imageRBG(i,j,1:3);
        newImg(i,j,1:3) = [0,128,128]'+Tycbcr*n;
    end
end

while mod(size(newImg,1),16)~=0
    newImg(end,:,:)= [];
end
while mod(size(newImg,2),16)~=0
    newImg(:,end,:)= [];
end

imageY  = newImg(:,:,1);

if isequal(subimg,[4,4,4])
    
    imageCb = newImg(:,:,2);
    imageCr = newImg(:,:,3);
    
elseif isequal(subimg,[4,2,2])
    
    imageCb = newImg(:,1:2:end,2);
    imageCr = newImg(:,1:2:end,3);
    
elseif isequal(subimg,[4,2,0])
    
    imageCb = newImg(1:2:end,1:2:end,2);
    imageCr = newImg(1:2:end,1:2:end,3);
    
else
    fprintf("Subsampling Vector: ")
    disp(subimg)
    error("Invalid Subsampling Vector. Use [4 4 4] or [4 2 2] or [4 2 0] instead.")
end
end

%% Convert YCbCr Image to RGB Image
function imageRGB = convert2rgb(imageY, imageCb, imageCr, subimg)

% Preproccesing
imageYCbCr(:,:,1)=imageY;

if isequal(subimg,[4,4,4])
    
    imageYCbCr(:,:,2) = imageCb;
    imageYCbCr(:,:,3) = imageCr;
    
elseif isequal(subimg,[4,2,2])
    
    imageYCbCr(:,:,2) = kron(imageCb,ones(1,2));
    imageYCbCr(:,:,3) = kron(imageCr,ones(1,2));
    
elseif isequal(subimg,[4,2,0])
    
    imageYCbCr(:,:,2) = kron(imageCb,ones(2));
    imageYCbCr(:,:,3) = kron(imageCr,ones(2));
    
else
    fprintf("Subsampling Vector: ")
    disp(subimg)
    error("Invalid Subsampling Vector. Use [4 4 4] or [4 2 2] or [4 2 0] instead.")
end

% Transformation Matrix
Trgb = [1,0,1.402;
    1,-.344136,-.714136;
    1,1.772,0];

imageRGB = zeros(size(imageYCbCr));
n = zeros(3,1);

% Transformation from YCbCr to RGB
for i = 1:size(imageYCbCr,1)
    for j = 1:size(imageYCbCr,2)
        n(:,1) = imageYCbCr(i,j,1:3);
        imageRGB(i,j,1:3) = Trgb*(n-[0,128,128]');
    end
end

imageRGB = uint8(imageRGB);

end

%% Discrete Cosine Transform of a block
function dctBlock = blockDCT(block)

dctBlock =  dct2(block);

end

%% Inverse Discrete Cosine Transform of a dct Block
function block = iBlockDCT(dctBlock)

block =  idct2(dctBlock);

end

%% Quantizer
function qBlock = quantizeJPEG(dctBlock, qTable, qScale)

qBlock = round(dctBlock./(qScale*qTable));

end

%% De-Quantizer
function dctBlock = dequantizeJPEG(qBlock, qTable, qScale)

dctBlock = qBlock.*(qScale*qTable);

end

%% Run Length Encoding
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

%% Inverse Run Length Encoding
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

%% JPEG Encoding
function JPEGenc = JPEGencode(img, subimg, qScale)

global Nx;
Tables;

YCbCr = cell(1,3);
[YCbCr{1},YCbCr{2},YCbCr{3}] = convert2ycbcr(img,subimg);

dim = size(YCbCr{1},1)/8+2*size(YCbCr{2},1)/8 + ...
    size(YCbCr{1},2)/8+2*size(YCbCr{2},2)/8;

JPEGenc = cell(dim+1,1);

tables.qTableL = qTable{1};
tables.qTableC = qTable{2};
tables.DCL = DCCategoryCode{1};
tables.DCC = DCCategoryCode{2};
tables.ACL = ACCategoryCode{1};
tables.ACC = ACCategoryCode{2};

JPEGenc{1} = tables;

count = 2;

for blockType = 1 : 3
    
    DCpred = 0;
    
    for indHor = 1 : size(YCbCr{blockType},1)/8
        for indVer = 1 : size(YCbCr{blockType},2)/8
            
            Nx.blkType = blockType;
            Nx.indHor = indHor;
            Nx.indVer = indVer;
            Nx.DCTable = DCCategoryCode{blockType};
            Nx.ACTable = ACCategoryCode{blockType};
            
            if blockType == 1
                qTable = tables.qTableL;
            else
                qTable = tables.qTableC;
            end
            
            block = YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8);
            dctBlock = blockDCT(block);
            
            qBlock = quantizeJPEG(dctBlock,qTable,qScale);
            
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
            
            Nx.huffStream = huffEnc(runSymbols);
            
            JPEGenc{count} = Nx;
            count = count + 1;
        end
    end
end

end

%% JPEG Decode
function imgRec = JPEGdecode(JPEGenc, subimg, qScale)

global DCTable
global ACTable

YCbCr = cell(1,3);

for count = 2 : length(JPEGenc)
    
    blockType = JPEGenc{count}.blkType;
    indHor  = JPEGenc{count}.indHor;
    indVer  = JPEGenc{count}.indVer;
    
    DCTable = JPEGenc{count}.DCTable;
    ACTable = JPEGenc{count}.ACTable;
     
    if blockType == 1
                qTable = JPEGenc{1}.qTableL;
            else
                qTable = JPEGenc{1}.qTableC;
     end
            
    if indHor == 1 && indVer == 1
        DCpred = 0;
    end
    
    runSymbols = huffDec(JPEGenc{count}.huffStream);
    
    qBlock = irunLength(runSymbols,DCpred);
    
    DCpred = qBlock(1,1);
    
    dctBlock = dequantizeJPEG(qBlock,qTable,qScale);
    block = iBlockDCT(dctBlock);
    
    YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8) = block;
    
end

imgRec = convert2rgb(YCbCr{1},YCbCr{2},YCbCr{3}, subimg);

end