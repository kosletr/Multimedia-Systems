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
% i=130; j=130;
Block = originalImage(i+1:i+8,j+1:j+8,1);
dctBlock = blockDCT(Block);
%%%
qBlock = quantizeJPEG(dctBlock,qTableY,1);
% qBlock = dctBlock;
%%%
runSymbols = runlength(qBlock,0);


qblock = irunlength(runSymbols,0);
%%%
dctblock = dequantizeJPEG(qBlock,qTableY,1);
% dctblock = qblock;
%%%
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

function JPEGenc = JPEGencode(img,subimg,qScale)



end