% Konstantinos Letros 8851
% Multimedia Systems Project
% JPEG Encoding

function JPEGenc = JPEGencode(img, subimg, qScale)

global blockType;

global DCCategoryCode;
global ACCategoryCode;
Tables;


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
                indVer = mod(k-1,2)+1 + mod(2*floor((size(YCbCr{blockType},2)/8*(j-1)+k-1)/4),size(YCbCr{blockType},2)/8);
%               fprintf('(%d,%d) --> (%d,%d) \n',j,k,indHor,indVer)
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