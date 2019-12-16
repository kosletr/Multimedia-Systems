% Konstantinos Letros 8851
% Multimedia Systems Project
% MSE Calculation

%% Clean the screen
close all
clc
clear
%% Main

img = imread('pics/flower2.jpg');
resultsFunc(img)

%% Results Function
function resultsFunc(img)

figure
imshow(img)

subimg = [4,4,4;
    4,2,2;
    4,2,0];

qScale = 1;

bitsNum = zeros(length(qScale),1);
mseImg = length(qScale);

for s = 1 : 3
    
    JPEGenc = JPEGencode(img,subimg(s,:),qScale);
    
    for c = 2 : length(JPEGenc)
        bitsNum(q) = bitsNum(q) + length(JPEGenc{c}.huffStream)*8;
    end
    
    imgREC = JPEGdecode(JPEGenc,subimg(s,:),qScale);
    %         figure
    %         imshow(imgREC)
    %
    %         imgRec = zeros(size(img));
    %         imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;
    %         figure
    %         imshow(uint8(double(img)-imgRec))
    
    mseImg(q) = mseEval(img,imgREC);
    
    
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

%% Save Plot
function savePlot(name)

% Resize current figure to fullscreen for higher resolution image
set(gcf, 'Position', get(0, 'Screensize'));

% Save current figure with the specified name
saveas(gcf, join([name,'.jpg']));

% Resize current figure back to normal
set(gcf,'position',get(0,'defaultfigureposition'));

end