% Konstantinos Letros 8851
% Multimedia Systems Project
% Demo 2

%% Clean The Screen
clc
clear
clear global
close all

%% Load Images

imgStruct = load('img1_down.mat');
img1 = imgStruct.img1_down;
imgStruct = load('img2_down.mat');
img2 = imgStruct.img2_down;

[rgbEntropy1,qBlockEntropy1,rleEntropy1] = EntropyMeasure(img1,[4 2 2],0.6)
[rgbEntropy2,qBlockEntropy2,rleEntropy2] = EntropyMeasure(img2,[4 4 4],5)

%% Function to measure Entropy
function [rgbEntropy,qDCTEntropy,rleEntropy] = EntropyMeasure(img,subimg,qScale)

rgbEntropy = myEntropy(img(:,:,1))*numel(img(:,:,1))+ ...
    myEntropy(img(:,:,2))*numel(img(:,:,2))+myEntropy(img(:,:,3))*numel(img(:,:,3));

qDCTEntropy = 0;
rleEntropy = 0;

YCbCr = cell(1,3);
[YCbCr{1},YCbCr{2},YCbCr{3}] = convert2ycbcr(img,subimg);
global zigZagTable;
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
    
    rleSymbols = rleStack(:,1)*1e4+rleStack(:,2);
    
    qDCTEntropy = qDCTEntropy + myEntropy(qDCTStack)*numel(qDCTStack);
    rleEntropy = rleEntropy + myEntropy(rleSymbols)*numel(rleSymbols);
    
end

end

function res = myEntropy(Array)

tab = tabulate(Array(:));
p = tab(:,end)/100;
en = p.*log2(1./p);
en(isnan(en)) = 0;
res = sum(en);

end