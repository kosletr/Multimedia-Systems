% Konstantinos Letros 8851
% Multimedia Systems Project
% JPEG Encoding

function JPEGenc = JPEGencode(img, subimg, qScale)

global Nx;
Tables;

YCbCr = cell(1,3);
[YCbCr{1},YCbCr{2},YCbCr{3}] = convert2ycbcr(img,subimg);

dim = size(YCbCr{1},1)/8+2*size(YCbCr{2},1)/8 + ...
    size(YCbCr{1},2)/8+2*size(YCbCr{2},2)/8;

JPEGenc = cell(dim+1,1);

tables.qTableL = qTable{1};
tables.qTableC = qTable{2};
tables.DCL = DCCategoryCode{1};
tables.DCC = DCCategoryCode{2};
tables.ACL = ACCategoryCode{1};
tables.ACC = ACCategoryCode{2};

JPEGenc{1} = tables;

count = 2;

for blockType = 1 : 3
    
    DCpred = 0;
    
    for indHor = 1 : size(YCbCr{blockType},1)/8
        for indVer = 1 : size(YCbCr{blockType},2)/8
            
            Nx.blkType = blockType;
            Nx.indHor = indHor;
            Nx.indVer = indVer;
            Nx.DCTable = DCCategoryCode{blockType};
            Nx.ACTable = ACCategoryCode{blockType};
            
            if blockType == 1
                qTable = tables.qTableL;
            else
                qTable = tables.qTableC;
            end
            
            block = YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8);
            dctBlock = blockDCT(block);
            
            qBlock = quantizeJPEG(dctBlock,qTable,qScale);
            
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
            
            Nx.huffStream = huffEnc(runSymbols);
            
            JPEGenc{count} = Nx;
            count = count + 1;
        end
    end
end

end