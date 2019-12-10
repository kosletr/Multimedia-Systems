% Konstantinos Letros 8851
% Multimedia Systems Project

%% Clean the screen

clc
clear
close all;

%% Defining Tables

global qTableY
global qTableCbCr
global zigZagTable
Tables

%% Testing

originalImage = imread('flower.jpg');
% figure
% imshow(originalImage)
%
% subVec = [4,2,0];
% [Y,Cb,Cr] = convert2ycbcr(originalImage,subVec);
% RGBimage = convert2RGB(Y,Cb,Cr, subVec);
%
% figure
% imshow(uint8(RGBimage))

i=randi([0,239]);
j=randi([0,134]);
Block = originalImage(i+1:i+8,j+1:j+8,1);
dctBlock = blockDCT(Block);
qBlock = quantizeJPEG(dctBlock,qTableY,1);
runSymbols = runlength(qBlock,0);
huffstream = huffEnc(runSymbols,1);

runSymbols = huffDec(huffstream,1);
qblock = irunlength(runSymbols,0);
dctblock = dequantizeJPEG(qBlock,qTableY,1);
block = iblockDCT(dctblock);
block = uint8(block);

% sum(sum(qblock==qBlock))/64

%% Convert RGB Image to YCbCr Image
function [ImageY,ImageCb,ImageCr] = convert2ycbcr(imageRBG,subimg)

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

ImageY  = newImg(:,:,1);

if subimg == [4,4,4]
    
    ImageCb = newImg(:,:,2);
    ImageCr = newImg(:,:,3);
    
elseif subimg == [4,2,2]
    
    ImageCb = newImg(:,1:2:end,2);
    ImageCr = newImg(:,1:2:end,3);
    
elseif subimg == [4,2,0]
    
    ImageCb = newImg(1:2:end,1:2:end,2);
    ImageCr = newImg(1:2:end,1:2:end,3);
    
else
    fprintf("Subsampling Vector: ")
    disp(subimg)
    error("Invalid Subsampling Vector. Use [4 4 4] or [4 2 2] or [4 2 0] instead.")
end
end

%% Convert YCbCr Image to RGB Image
function ImageRGB = convert2RGB(ImageY,ImageCb, ImageCr, subimg)

% Preproccesing
imageYCbCr(:,:,1)=ImageY;

if subimg == [4,4,4]
    
    imageYCbCr(:,:,2) = ImageCb;
    imageYCbCr(:,:,3) = ImageCr;
    
elseif subimg == [4,2,2]
    
    imageYCbCr(:,:,2) = kron(ImageCb,ones(1,2));
    imageYCbCr(:,:,3) = kron(ImageCr,ones(1,2));
    
elseif subimg == [4,2,0]
    
    imageYCbCr(:,:,2) = kron(ImageCb,ones(2));
    imageYCbCr(:,:,3) = kron(ImageCr,ones(2));
    
else
    fprintf("Subsampling Vector: ")
    disp(subimg)
    error("Invalid Subsampling Vector. Use [4 4 4] or [4 2 2] or [4 2 0] instead.")
end

% Transformation Matrix
Trgb = [1,0,1.402;
    1,-.344136,-.714136;
    1,1.772,0];

ImageRGB = zeros(size(imageYCbCr));
n = zeros(3,1);

% Transformation from YCbCr to RGB
for i = 1:size(imageYCbCr,1)
    for j = 1:size(imageYCbCr,2)
        n(:,1) = imageYCbCr(i,j,1:3);
        ImageRGB(i,j,1:3) = Trgb*(n-[0,128,128]');
    end
end


end

%% Discrete Cosine Transform of a block
function dctBlock = blockDCT(block)

dctBlock =  dct2(block);

end

%% Inverse Discrete Cosine Transform of a dct Block
function block = iblockDCT(dctBlock)

block =  idct2(dctBlock);

end

%% Quantizer
function qBlock = quantizeJPEG(dctBlock,qTable,qScale)

qBlock = round(dctBlock./(qScale*qTable));

end

%% De-Quantizer
function dctBlock = dequantizeJPEG(qBlock,qTable,qScale)

dctBlock = qBlock.*(qScale*qTable);

end

%% Run Length Encoding
function runSymbols = runlength(qBlock,DCpred)

% Part 1
global zigZagTable
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
i = 1;
k = 1;

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

function qBlock = irunlength(runSymbols, DCpred)

% Part 1
runSymbols(end,2) = NaN; % Adding ending delimeter
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
qBlock = zeros(sqrt(length(vec)));

% Place values back, in zig-zag form
for i = 1:length(vec)
    [r,c]=find(zigZagTable==i);
    qBlock(r,c) = vec(i);
end

qBlock(1,1) = qBlock(1,1) + DCpred;

end

function huffstream = huffEnc(runSymbols,y)

global DCCategoryCode
global ACCategoryCode

binaryStream = [];

for k = 1:length(runSymbols)

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
    DCvector=[DCCategoryCode{y}{huffmanInd}-'0',magnitude];
    binaryStream = [binaryStream,DCvector];
else
    % Find Index of Run/Category in AC table
    huffmanInd = 10*runSymbols(k,1)+(category+1);
    if runSymbols(k,1)==15 % index 152 --> ZRL
        huffmanInd = huffmanInd + 1;
    elseif runSymbols(k,1) == 0 && runSymbols(k,2) == 0 % EOB
        magnitude = [];
    end
    ACvector = [ACCategoryCode{y}{huffmanInd}-'0',magnitude];
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
huffstream = zeros(1,length(binaryStream)/8);
while i < length(binaryStream)
   huffstream(k) = uint8(bi2de(binaryStream(i+1:i+8),'left-msb'));
   i = i + 8;
   k = k + 1;
end

end

function runSymbols = huffDec(huffstream,y)

global DCCategoryCode
global ACCategoryCode

% Convert bytes to binary stream
binaryStream = [];
for i = 1:length(huffstream)
binaryStream = [binaryStream,de2bi(huffstream(i),8,'left-msb')];
end

% Initialize
k=1;
bit = 1;
stop = 2;

% Convert binary vector to string
strStream = sprintf('%d',binaryStream(bit:stop));

% DC Difference Coefficient
% Repeat until unique reference (DC Table)
while length(find(DCCategoryCode{y}==strStream))~=1
    stop = stop + 1;
    strStream = sprintf('%d',binaryStream(bit:stop));
end

category = find(DCCategoryCode{y}==strStream)-1;
bit = stop + 1;
stop = stop + category;

% 0 -> Negative Sign
if binaryStream(bit)==0
    DCdiff = bi2de(~binaryStream(bit:stop),'left-msb');
else
    DCdiff = bi2de(binaryStream(bit:stop),'left-msb');
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
    while length(find(ACCategoryCode{y}==strStream))~=1
        stop = stop + 1;
        strStream = sprintf('%d',binaryStream(bit:stop));
    end
    
    if isequal(binaryStream(bit:stop) , [1 0 1 0]) % EOB
        runSymbols(k,:)= [0,0];
        break;
    end
    
    indexAC = find(ACCategoryCode{y}==strStream)-1;
    precZeros = fix(indexAC/10);
    category = mod(indexAC,10);
    if precZeros == 15
        category = category - 1;
    end
    
    bit = stop + 1;
    stop = stop + category;
    
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

function JPEGenc = JPEGencode(img,subimg,qScale)



end