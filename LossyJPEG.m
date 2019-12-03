% Konstantinos Letros 8851
% Multimedia Systems Project

%% Clean the screen

clc
clear
close all;

%% Testing

originalImage = imread('flower.jpg');
figure
imshow(originalImage)

subVec = [4,2,0];
[Y,Cb,Cr] = convert2ycbcr(originalImage,subVec);
RGBimage = convert2RGB(Y,Cb,Cr, subVec);

figure
imshow(uint8(RGBimage))

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

qBlock = round(dctBlock./(qScale*qTable));

end

%% De-Quantizer
function dctBlock = dequantizeJPEG(qBlock,qTable,qScale)

dctBlock = qBlock.*(qScale*qTable);

end

%% Zig Zag Array
function Z = zigzag(N)

Z = zeros(N);

% First Row
for j = 2:N
    if mod((1+j),2)==0
        Z(1,j) = Z(1,j-1)+ 2*(j - 1);
    else
        Z(1,j) = Z(1,j-1) + 1;
    end
end

% First Element of the Last Row
if(mod(N,2)==0)
    Z(N,1) = Z(1,N)+(N-1);
else
    Z(N,1) = Z(1,N)-(N-1);
end

% Last Row
for j = 2:N
    if mod((N+j),2)==1
        Z(N,j) = Z(N,j-1)+ 2*(N+1-j);
    else
        Z(N,j) = Z(N,j-1) + 1;
    end
end

% Upper Left Array
for i = 2:N
    for j = 1:N-i
        if(mod((i+j),2)==1)
            Z(i,j)=Z(i-1,j+1)+1;
        else
            Z(i,j)=Z(i-1,j+1)-1;
        end
    end
end

% Lower Right Array
for i = N-1:-1:2
    for j = N-i+1:N
        if(mod((i+j),2)==1)
            Z(i,j)=Z(i+1,j-1)-1;
        else
            Z(i,j)=Z(i+1,j-1)+1;
        end
    end
end

end