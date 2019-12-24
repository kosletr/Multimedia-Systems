% Konstantinos Letros 8851
% Multimedia Systems Project
% MSE Calculation

%% Clean the screen
close all
clc
clear
%% Main

imgStruct = load('img1_down.mat');
img1 = imgStruct.img1_down;
imgStruct = load('img2_down.mat');
img2 = imgStruct.img2_down;

resultsFunc(img1)
resultsFunc(img2)

%% Results Function
function resultsFunc(img)

figure
imshow(img)

subimg = [4,4,4;
    4,2,2;
    4,2,0];

qScale = [0.1; 0.3; 0.6; 1; 2; 5; 10];

bitsNum = zeros(length(qScale),1);
mseImg = length(qScale);
imageSize = size(img,1)*size(img,2)*size(img,3)*8;

for s = 1 : 3
    for q = 1 :length(qScale)
        
        JPEGenc = JPEGencode(img,subimg(s,:),qScale(q));
        
        for c = 2 : length(JPEGenc)
            bitsNum(q) = bitsNum(q) + length(JPEGenc{c}.huffStream)*8;
        end
        
        
        compRatio(q) = bitsNum(q)/imageSize;
        
        imgREC = JPEGdecode(JPEGenc,subimg(s,:),qScale(q));
        %         figure
        %         imshow(imgREC)
        %
        %         imgRec = zeros(size(img));
        %         imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;
        %         figure
        %         imshow(uint8(double(img)-imgRec))
        
        mseImg(q) = mseEval(img,imgREC);
        
    end
    
    fprintf("Compression Ratios: %f", compRatio)
    
    figure
    bar(qScale,mseImg)
    title(['Subsampling: [',num2str(subimg(s,:)),']'])
    xlabel('qScale Values')
    ylabel('Mean Sqaure Error')
    
    figure
    plot(bitsNum,mseImg)
    title(['Subsampling: [',num2str(subimg(s,:)),']'])
    xlabel('Number of Bits in Bitstream')
    ylabel('Mean Sqaure Error')
end
end

%% Evaluate MSE
function MSE = mseEval(img,imgREC)

% Zero Padding to make compatible dimensions
imgRec = zeros(size(img));
imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;

% Calculate MSE
MSE  = 1/((size(img,1)*size(img,2)) * sum(sum(sum((double(img)-imgRec).^2))));

end