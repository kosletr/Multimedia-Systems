originalImage = imread('pics/flower2.jpg');
figure
imshow(originalImage)

subimg = [4,2,0];
qScale = 1;

JPEGenc = JPEGencode(originalImage,subimg,qScale);
imgREC = JPEGdecode(JPEGenc,subimg,qScale);
figure
imshow(imgREC)