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

[rgbEntropy,qBlockEntropy,rleEntropy] = EntropyMeasure(img1,[4 2 0],1)

%% Function to measure Entropy
function [rgbEntropy,qDCTEntropy,rleEntropy] = EntropyMeasure(img,subimg,qScale)

rgbEntropy = entropy(img(:,:,1))+entropy(img(:,:,2))+entropy(img(:,:,3));

qDCTEntropy = 0;
rleEntropy = 0;

YCbCr = cell(1,3);
[YCbCr{1},YCbCr{2},YCbCr{3}] = convert2ycbcr(img,subimg);
Tables;

for blockType = 1 : 3
    
    qDCTStack = [];
    rleStack = [];
    DCpred = 0;
    
    for indHor = 1 : size(YCbCr{blockType},1)/8
        for indVer = 1 : size(YCbCr{blockType},2)/8
            
            block = YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8);
            dctBlock = blockDCT(block);
            
            qBlock = quantizeJPEG(dctBlock,qTable{blockType},qScale);
            qDCTStack = [qDCTStack ; qBlock];
            
            
            runSymbols = runLength(qBlock,DCpred);
            rleStack = [rleStack ; runSymbols];
            DCpred = qBlock(1,1);
            
        end
    end
    qDCTEntropy = qDCTEntropy + entropy(qDCTStack);
    rleEntropy = rleEntropy + entropy(rleStack);
            
            figure
            histogram(qDCTStack)
            title(num2str(entropy(qDCTStack)))
            figure
            histogram(rleStack)
            title(num2str(entropy(rleStack)))
end

end