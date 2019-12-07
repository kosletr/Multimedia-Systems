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

% Quantization Tables

qTableY = [ ...
    16 11 10 16 124 140 151 161;
    12 12 14 19 126 158 160 155;
    14 13 16 24 140 157 169 156;
    14 17 22 29 151 187 180 162;
    18 22 37 56 168 109 103 177;
    24 35 55 64 181 104 113 192;
    49 64 78 87 103 121 120 101;
    72 92 95 98 112 100 103 199];

qTableCbCr = [ ...
    17 18 24 47 99 99 99 99;
    18 21 26 66 99 99 99 99;
    24 26 56 99 99 99 99 99;
    47 66 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99];

% ZigZag Table 8x8

zigZagTable = [ ...
    1  2  6  7  15 16 28 29;
    3  5  8  14 17 27 30 43;
    4  9  13 18 26 31 42 44;
    10 12 19 25 32 41 45 54;
    11 20 24 33 40 46 53 55;
    21 23 34 39 47 52 56 61;
    22 35 38 48 51 57 60 62;
    36 37 49 50 58 59 63 64];

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



qblock = irunlength(runSymbols,0);
dctblock = dequantizeJPEG(qBlock,qTableY,1);
block = iblockDCT(dctblock);
block = uint8(block);


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
zzBlockVec = [ NaN ; zzBlockVec ; NaN]; % Adding delimeters
i = 1;
k = 1;

% RLE

while i<length(zzBlockVec)
    i = i + 1;
    countZeros = 0;
    j = i;
    while zzBlockVec(j-1)==0
        countZeros = countZeros + 1;
        j = j -1;
    end
    if zzBlockVec(i)~=0
        runSymbols(k,1)=countZeros;
        runSymbols(k,2)=zzBlockVec(i);
        k = k + 1;
    end
end
runSymbols(end,2)=0; % Removing ending delimeter

% max length until reset = 15
i = 1;
while i <= length(runSymbols)
   if runSymbols(i,1)>15
        runSymbols = [runSymbols(1:i-1,:);[15,0];runSymbols(i,:)-[16,0];runSymbols(i+1:end,:)];
         i = 1;
   end
   i = i + 1;
end

end

function qBlock = irunlength(runSymbols, DCpred)

% Part 1
runSymbols(end,2) = NaN; % Adding ending delimeter
vec = [];

% Inverse RLE
for i=1:length(runSymbols)
    vec = [vec;zeros(runSymbols(i,1),1)];
    vec(end+1,1)=runSymbols(i,2);
end

vec(end)=[]; % Removing ending delimeter

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