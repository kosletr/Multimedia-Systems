S = dir('*.jpg');
for k = 1:numel(S)
    img = imread(S(k).name);
    imageCrop(img,S(k).name)
end



function imageCrop(img,filename)

col = [];
for j = 1:size(img,2)
    col(j)=mean(img(:,j,1));
    if col(j)~=255
        break
    end
end

left=j;
col = [];
for j = size(img,2):-1:1
    col(j)=mean(img(:,j,1));
    if col(j)~=255
        break
    end
end

right = j;

row = [];
for i = 1:size(img,1)
    row(i)=mean(img(i,:,1));
    if row(i)~=255
        break
    end
end

up = i;

row = [];
for i = size(img,1):-1:1
    row(i)=mean(img(i,:,1));
    if row(i)~=255
        break
    end
end

down = i;

imwrite(img(up:down,left:right,:),filename);
end