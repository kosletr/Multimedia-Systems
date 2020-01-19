% Konstantinos Letros 8851
% Multimedia Systems Project
% qTables Distortion

%% Clean the screen

close all
clc
clear
format long g

%% Main

% Load Images
imgStruct = load('img1_down.mat');
img1 = imgStruct.img1_down;
imgStruct = load('img2_down.mat');
img2 = imgStruct.img2_down;

% Results Function
% resultsBFunc(img1,[4 2 2],1)
resultsBFunc(img2,[4 4 4],1)

% Save Plots
h =  findobj('type','figure');
for i = 1 : length(h)
    figure(i)
    savePlot([mfilename,'_',num2str(i)])
end

%% Results Function
function resultsBFunc(img,subimg,qScale)

% Show Original Image
figure
imshow(img)
title('Original Image')

terms = [0;20;40;50;60;63];

% Initialization
bitsNum = zeros(length(terms),1);
mseImg = zeros(length(terms),1);
compRatio = zeros(length(terms),1);

% Original Image Size
imageSize = size(img,1)*size(img,2)*size(img,3);

fprintf("Image Results \n\n")

% Main Loop
for n = 1 : length(terms)
    
    JPEGenc = JPEGencodeD(img,subimg,qScale,n);
    
    % Count bits of bit-Stream
    for c = 2 : length(JPEGenc)
        bitsNum(n) = bitsNum(n) + length(JPEGenc{c}.huffStream)*8;
    end
    
    % Calculate Compression Ratio
    compRatio(n) = (imageSize*8)/bitsNum(n);
    
    imgREC = JPEGdecode(JPEGenc);
    
    % Show reconstructed Image
    figure
    imshow(imgREC)
    title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale),...
        ' - AC Terms Removed: ',num2str(terms(n)),' - Image Reconstruction'])
    
    
    % Zero Padding to make compatible dimensions
    imgRec = zeros(size(img));
    imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;
    
    % Show Error Image
    figure
    imshow(uint8(double(img)-imgRec))
    title(['AC Terms Removed: ',num2str(terms(n)),' - Subsampling: [',...
        num2str(subimg),']',' - qScale: ',num2str(qScale),' - Image Error'])
    
    % Evaluate Mean Square Error
    mseImg(n) = 1/imageSize * sum((double(img(:))-imgRec(:)).^2);
    
    fprintf("AC Terms Removed: %d - Compression Ratio: %f - MSE: %f - Number of bits: %d \n", ...
        terms(n), compRatio(n), mseImg(n), bitsNum(n))
    
end

fprintf("\n\n")


figure
bar(terms,mseImg)
title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale)])
xlabel('Number of Deleted qTable Terms')
ylabel('Mean Sqaure Error')

figure
plot(bitsNum,mseImg,'x-','MarkerSize',10,'LineWidth',1.5)
title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale)])
xlabel('Number of Bits in Bitstream')
ylabel('Mean Sqaure Error')

end

%% JPEG Encode Distorted qTables


function JPEGenc = JPEGencodeD(img, subimg, qScale, n)

global blockType;

global DCCategoryCode;
global ACCategoryCode;
Tables;
qTable = distortedqTables(n);


types = ["Y","Cb","Cr"];

YCbCr = cell(1,3);
[YCbCr{1},YCbCr{2},YCbCr{3}] = convert2ycbcr(img,subimg);

dim = size(YCbCr{1},1)/8+2*size(YCbCr{2},1)/8 + ...
    size(YCbCr{1},2)/8+2*size(YCbCr{2},2)/8;

JPEGenc = cell(dim+1,1);

tables.qTableL = qScale*qTable{1};
tables.qTableC = qScale*qTable{2};
tables.DCL = DCCategoryCode{1};
tables.DCC = DCCategoryCode{2};
tables.ACL = ACCategoryCode{1};
tables.ACC = ACCategoryCode{2};

JPEGenc{1} = tables;

count = 2;

for blockType = 1 : 3
    
    DCpred = 0;
    
    for j = 1 : size(YCbCr{blockType},1)/8
        for k = 1 : size(YCbCr{blockType},2)/8
            
            if isequal(subimg,[4 2 0]) && blockType == 1
                indHor = mod(ceil((size(YCbCr{blockType},2)/8*(j-1)+k)/2)+1,2)+1+2*floor((j-1)/2);
                indVer = mod(k-1,2)+1 + ...
                    mod(2*floor((size(YCbCr{blockType},2)/8*(j-1)+k-1)/4),size(YCbCr{blockType},2)/8);
            else
                indHor = j;
                indVer = k;
            end
            
            
            blk.blkType = types(blockType);
            blk.indHor = indHor;
            blk.indVer = indVer;
            
            block = YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8);
            dctBlock = blockDCT(block);
            
            qBlock = quantizeJPEG(dctBlock,qTable{blockType},qScale);
            
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
            
            blk.huffStream = huffEnc(runSymbols);
            
            JPEGenc{count} = blk;
            count = count + 1;
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