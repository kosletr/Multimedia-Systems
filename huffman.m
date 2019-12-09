
h=huffEnc(runSymbols,1)
runSymbols
length(h)

function huffstream = huffEnc(runSymbols,y)

global DCCategoryCode
global ACCategoryCode
Tables

huffstream = [];

for k = 1:length(runSymbols)

% Find the Category of the non-zero DCT coeffs
category = 0;

if runSymbols(k,2)~=0
    while category<11
        if 2^(category-1)<=abs(runSymbols(k,2)) && abs(runSymbols(k,2))<=(2^category)-1
            break
        end
        category = category + 1;
    end
else
    category = 0;
end

% Convert non-zero DCT Coeffs to binary

if runSymbols(k,2)<0 % 1's Complement
   magnitude = double(~de2bi(-runSymbols(k,2),'left-msb')); 
else
    magnitude = de2bi(runSymbols(k,2),'left-msb');
end


if k==1
    % Find Index of Category in DC table
    huffmanInd = category + 1;
    DCvector=[DCCategoryCode{y}{huffmanInd}-'0',magnitude]
    huffstream = [huffstream,DCvector];
else
    % Find Index of Run/Category in AC table
    huffmanInd = 10*runSymbols(k,1)+(category+1);
    if runSymbols(k,1)==15 % index 152 --> ZRL
        huffmanInd = huffmanInd + 1;
    elseif runSymbols(k,1) == 0 && runSymbols(k,2) == 0 % EOB
        magnitude = [];
    end
    ACvector = [ACCategoryCode{y}{huffmanInd}-'0',magnitude]
    huffstream = [huffstream,ACvector];
end

end

end