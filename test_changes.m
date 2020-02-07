clear all;
close all;
clc

load('img1_down.mat')

load('img2_down.mat')

test1 = imread('Pics/test1.jpeg');

index = imread('Pics/index.jfif');

test2 = imread('Pics/test2.jfif');

test3 = imread('Pics/test3.jfif');

test4 = imread('Pics/wUn6JO1.jpg');

cellA = cell(7,1);
cellA{1,1} = img1_down; cellA{2,1} = img2_down; cellA{3,1} = test1; cellA{4,1} = test2;
cellA{5,1} = test3; cellA{6,1} = test4; cellA{7,1} = index;
for i = 1:length(cellA)
    JPEGenc = JPEGencode(cellA{i,1}, [4 4 4], 1);
    imgRec = JPEGdecode(JPEGenc);
    figure();
    subplot(1,2,1)
    imshow(cellA{i,1})
    title('Original Image', 'Interpreter', 'latex')
    subplot(1,2,2)
    imshow(imgRec);
    title_str = ['Reconstructed Image'];
    title(title_str, 'Interpreter', 'latex')
    pause(0.1);
end