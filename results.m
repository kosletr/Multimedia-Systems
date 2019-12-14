% Konstantinos Letros 8851
% Multimedia Systems Project
% MSE Calculate

originalImage = imread('pics/flower2.jpg');
figure
imshow(originalImage)

subimg = [4,2,0];
qScale = 1;

JPEGenc = JPEGencode(originalImage,subimg,qScale);
imgREC = JPEGdecode(JPEGenc,subimg,qScale);
figure
imshow(imgREC)

mseEval(originalImage,imgREC)

function MSE = mseEval(img,imgREC)

% Zero Padding to make compatible dimensions
imgRec = zeros(size(img));
imgRec(1:size(imgREC,1),1:size(imgREC,2),1:size(imgREC,3)) = imgREC;

% Calculate MSE
MSE  = 1/((size(img,1)*size(img,2)) * sum(sum(sum((double(img)-imgRec).^2))));

end