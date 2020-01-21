% Konstantinos Letros 8851
% Multimedia Systems Project
% Convert RGB Image to YCbCr Image

function [imageY, imageCb, imageCr] = convert2ycbcr(imageRGB, subimg)

% Transformation Matrix
Tycbcr = [.299,.587,.114;
    -.168736,-.331264,.5;
    .5,-.418688,-.081312];

newImg = zeros(size(imageRGB));
n = zeros(3,1);

% Transformation from RGB to YCbCr
for i = 1:size(imageRGB,1)
    for j = 1:size(imageRGB,2)
        n(:,1) = imageRGB(i,j,1:3);
        newImg(i,j,1:3) = [0,128,128]'+Tycbcr*n;
    end
end

% Remove rows, cols to fit 8x8 blocks
while mod(size(newImg,1),16)~=0
    newImg(end,:,:)= [];
end
while mod(size(newImg,2),16)~=0
    newImg(:,end,:)= [];
end

% Sub-sampling
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