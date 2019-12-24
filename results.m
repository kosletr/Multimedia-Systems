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

resultsFunc(img1,[4 2 2])
resultsFunc(img2,[4 4 4])

h =  findobj('type','figure');
for i = 1 : length(h)
    figure(i)
    savePlot([mfilename,'_',num2str(i)])
end

%% Results Function
function resultsFunc(img,subimg)

figure
imshow(img)

qScale = [0.1; 0.3; 0.6; 1; 2; 5; 10];

bitsNum = zeros(length(qScale),1);
mseImg = zeros(length(qScale),1);
imageSize = size(img,1)*size(img,2)*size(img,3)*8;
compRatio = zeros(length(qScale),1);

fprintf("\n\n Image Results \n\n")

for q = 1 :length(qScale)
    
    JPEGenc = JPEGencode(img,subimg,qScale(q));
    
    for c = 2 : length(JPEGenc)
        bitsNum(q) = bitsNum(q) + length(JPEGenc{c}.huffStream)*8;
    end
    
    
    compRatio(q) = imageSize/bitsNum(q);
    
    imgREC = JPEGdecode(JPEGenc,subimg,qScale(q));
    figure
    imshow(imgREC)
    title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale(q)),' - Image Reconstruction'])
    
    imgRec = zeros(size(img));
    imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;
    figure
    imshow(uint8(double(img)-imgRec))
    title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale(q)),' - Image Error'])
    
    
    mseImg(q) = mseEval(img,imgREC);
    
end

fprintf("qScale: %f - Compression Ratio: %f\n", qScale(q), compRatio)

figure
bar(qScale,mseImg)
title(['Subsampling: [',num2str(subimg),']'])
xlabel('qScale Values')
ylabel('Mean Sqaure Error')

figure
plot(bitsNum,mseImg)
title(['Subsampling: [',num2str(subimg),']'])
xlabel('Number of Bits in Bitstream')
ylabel('Mean Sqaure Error')

end

%% Evaluate MSE
function MSE = mseEval(img,imgREC)

% Zero Padding to make compatible dimensions
imgRec = zeros(size(img));
imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;

% Calculate MSE
MSE  = 1/((size(img,1)*size(img,2)) * sum(sum(sum((double(img)-imgRec).^2))));

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