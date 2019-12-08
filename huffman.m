function twoBytes = huffman(s)
global DCCategoryCode
global ACCategoryCode
Tables

%%
k=1;
% s = [1,-21];

category = 0;

if s(k,2)~=0
    while category<11
        if 2^(category-1)<=abs(s(k,2)) && abs(s(k,2))<=(2^category)-1
            break
        end
        category = category + 1;
    end
else
    category = 0;
end

num = zeros(1,category);

if s(k,2)<0
   s(k,2) = -s(k,2);
end
temp = dec2bin(s(k,2));
num(category-length(temp)+1:category) = temp - '0';

% DCvector=[DCCategoryCode{1}{category}-'0',num];
% twoBytes=uint16(bi2de(DCvector));

ACvector = [ACCategoryCode{1}{1+10*s(k,1)+category}-'0',num]
twoBytes=uint16(bi2de(ACvector));