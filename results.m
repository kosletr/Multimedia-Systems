% Konstantinos Letros 8851
% Multimedia Systems Project
% MSE Calculation

%% Clean the screen

close all
clc
clear
clear global
format long g

%% Main

% Load Images
imgStruct = load('img1_down.mat');
img1 = imgStruct.img1_down;
imgStruct = load('img2_down.mat');
img2 = imgStruct.img2_down;

% Results Function
resultsFunc(img1,[4 2 2])
resultsFunc(img2,[4 4 4])

% Save Plots

% h =  findobj('type','figure');
% for i = 1 : length(h)
%     figure(i)
%     savePlot([mfilename,'_',num2str(i)])
% end

%% Results Function
function resultsFunc(img,subimg)

N = mod(size(img,1),16);
M = mod(size(img,2),16);
img = img(1:end-N,1:end-M,:);

% Show Original Image
figure
imshow(img)
title('Original Image')


% qScale Values to be tested
qScale = [0.1; 0.3; 0.6; 1; 2; 5; 10];

% Initialization
bitsNum = zeros(length(qScale),1);
mseImg = zeros(length(qScale),1);
compRatio = zeros(length(qScale),1);

fprintf("Image Results \n\n")

% Main Loop
for q = 1 :length(qScale)
    
    JPEGenc = JPEGencode(img,subimg,qScale(q));
    
    % Count bits of bit-Stream
    for c = 2 : length(JPEGenc)
        bitsNum(q) = bitsNum(q) + length(JPEGenc{c}.huffStream);
    end
    
    % Calculate Compression Ratio
    compRatio(q) = (numel(img)*8)/bitsNum(q);
    
    imgREC = JPEGdecode(JPEGenc);
    
    % Show reconstructed Image
    figure
    imshow(imgREC)
    title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale(q)),' - Image Reconstruction'])
    
    
    % Zero Padding to make compatible dimensions
    imgRec = zeros(size(img));
    imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;
    
    % Show Error Image
    figure
    imshow(uint8(double(img)-imgRec))
    title(['Subsampling: [',num2str(subimg),']',' - qScale: ',num2str(qScale(q)),' - Image Error'])
    
    % Evaluate Mean Square Error
    mseImg(q) = 1/numel(img(:,:,1)) * sum(sum((double(img(:,:,1))-imgRec(:,:,1)).^2)) + ...
        1/numel(img(:,:,2)) * sum(sum((double(img(:,:,2))-imgRec(:,:,2)).^2)) + ...
        1/numel(img(:,:,3)) * sum(sum((double(img(:,:,3))-imgRec(:,:,3)).^2));
    
    fprintf("qScale: %f - Compression Ratio: %f - MSE: %f - Number of bits: %d \n", ...
        qScale(q), compRatio(q), mseImg(q), bitsNum(q))
    
end

fprintf("\n\n")


figure
bar(qScale,mseImg)
title(['Subsampling: [',num2str(subimg),']'])
xlabel('qScale Values')
ylabel('Mean Sqaure Error')

figure
plot(bitsNum,mseImg,'x-','MarkerSize',10,'LineWidth',1.5)
title(['Subsampling: [',num2str(subimg),']'])
xlabel('Number of Bits in Bitstream')
ylabel('Mean Sqaure Error')

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