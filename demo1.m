% Konstantinos Letros 8851
% Multimedia Systems Project
% Demo 1

%% Clean The Screen
clc
clear
close all

%% Load Images

imgStruct = load('img1_down.mat');
img1 = imgStruct.img1_down;
imgStruct = load('img2_down.mat');
img2 = imgStruct.img2_down;

%% Part A
part = 'A';

% Image 1
mainDemo(img1, [4 2 2], 1, part)

% Image 2
mainDemo(img2, [4 4 4], 2, part)

%% Part B
part = 'B';

% Image 1
mainDemo(img1, [4 2 2], 1, part)

% Image 2
mainDemo(img2, [4 4 4], 2, part)

%% Save Plots

h =  findobj('type','figure');
for i = 1 : length(h)
    figure(i)
    savePlot([mfilename,'_',num2str(i)])
end


%% Main Function

function mainDemo(img,subimg,imgNum,part)

[Y,Cb,Cr] = convert2ycbcr(img,subimg);

if part == 'B'
    
    qScale = [0.6; 5];
    
    YCbCr{1}=Y;
    YCbCr{2}=Cb;
    YCbCr{3}=Cr;
    
    partBDemo(YCbCr,qScale(imgNum));
end

RGBimage = convert2rgb(Y,Cr,Cb, subimg);

figure
imshow(img);
title('Original Image');

figure
imshow(RGBimage);
if part == "B"
    title(['Reconstructed Image - Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale(imgNum))])
else
    title('Reconstructed Image')
end

end

%% Part B Function

function YCbCrInv = partBDemo(YCbCr,qScale)
Tables;

for blockType = 1 : 3
    for indHor = 1 : size(YCbCr{blockType},1)/8
        for indVer = 1 : size(YCbCr{blockType},2)/8
            
            % Functions
            block = YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8);
            dctBlock = blockDCT(block);
            qBlock = quantizeJPEG(dctBlock,qTable{blockType},qScale);
            
            % Inverse Functions
            dctBlockInv = dequantizeJPEG(qBlock,qTable{blockType},qScale);
            blockInv = iBlockDCT(dctBlockInv);
            
            YCbCrInv{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8) = blockInv;
        end
    end
end

end

%% Function to automatically save plots in high resolution
function savePlot(name)

% Resize current figure to fullscreen for higher resolution image
set(gcf, 'Position', get(0, 'Screensize'));

% Save current figure with the specified name
saveas(gcf, join([name,'.jpg']));

% Resize current figure back to normal
set(gcf,'position',get(0,'defaultfigureposition'));

end