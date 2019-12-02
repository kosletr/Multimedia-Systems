% Konstantinos Letros 8851
% Multimedia Systems Project

%% Clean the screen

clc
clear
close all;

%% Testing

img = imread('flower2.jpg');
figure
imshow(img)

subVec = [4,2,0];
[Y,Cb,Cr] = convert2ycbcr(img,subVec);
img2 = convert2RGB(Y,Cb,Cr, subVec);

figure
imshow(uint8(img2))


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
    
    imageYCbCr(:,:,2) = kron(ImageCb,[1,1]);
    imageYCbCr(:,:,3) = kron(ImageCr,[1,1]);
    
elseif subimg == [4,2,0]
    
    imageYCbCr(:,:,2) = kron(ImageCb,[1,1;1,1]);
    imageYCbCr(:,:,3) = kron(ImageCr,[1,1;1,1]);
    
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

qBlock = 0;

end