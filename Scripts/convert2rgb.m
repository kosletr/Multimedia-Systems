%% Convert YCbCr Image to RGB Image
function imageRGB = convert2rgb(imageY, imageCb, imageCr, subimg)

% Preproccesing
imageYCbCr(:,:,1)=imageY;

if isequal(subimg,[4,4,4])
    
    imageYCbCr(:,:,2) = imageCb;
    imageYCbCr(:,:,3) = imageCr;
    
elseif isequal(subimg,[4,2,2])
    
    imageYCbCr(:,:,2) = kron(imageCb,ones(1,2));
    imageYCbCr(:,:,3) = kron(imageCr,ones(1,2));
    
elseif isequal(subimg,[4,2,0])
    
    imageYCbCr(:,:,2) = kron(imageCb,ones(2));
    imageYCbCr(:,:,3) = kron(imageCr,ones(2));
    
else
    fprintf("Subsampling Vector: ")
    disp(subimg)
    error("Invalid Subsampling Vector. Use [4 4 4] or [4 2 2] or [4 2 0] instead.")
end

% Transformation Matrix
Trgb = [1,0,1.402;
    1,-.344136,-.714136;
    1,1.772,0];

imageRGB = zeros(size(imageYCbCr));
n = zeros(3,1);

% Transformation from YCbCr to RGB
for i = 1:size(imageYCbCr,1)
    for j = 1:size(imageYCbCr,2)
        n(:,1) = imageYCbCr(i,j,1:3);
        imageRGB(i,j,1:3) = Trgb*(n-[0,128,128]');
    end
end

imageRGB = uint8(imageRGB);

end