% Konstantinos Letros 8851
% Multimedia Systems Project
% JPEG Decode

function imgRec = JPEGdecode(JPEGenc)
global blockType;
YCbCr = cell(1,3);
types = ["Y","Cb","Cr"];

for count = 2 : length(JPEGenc)
    
    blkType = JPEGenc{count}.blkType;
    indHor  = JPEGenc{count}.indHor;
    indVer  = JPEGenc{count}.indVer;
    
    blockType = find(blkType==types);
    
    if blockType == 1
        qTable = JPEGenc{1}.qTableL;
    else
        qTable = JPEGenc{1}.qTableC;
    end
    
    if indHor == 1 && indVer == 1
        DCpred = 0;
    end
    
    runSymbols = huffDec(double(JPEGenc{count}.huffStream));
    
    qBlock = irunLength(runSymbols,DCpred);
    
    DCpred = qBlock(1,1);
    
    dctBlock = dequantizeJPEG(qBlock,qTable,1);
    block = iBlockDCT(dctBlock);
    
    YCbCr{blockType}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8) = block;
    
end

lenYHor = size(YCbCr{1},1);
lenYVer = size(YCbCr{1},2);
lenCbCrHor = size(YCbCr{2},1);
lenCbCrVer = size(YCbCr{2},2);

if lenYHor/lenCbCrHor == 1 && lenYVer/lenCbCrVer == 1
    subimg = [4 4 4];
elseif lenYHor/lenCbCrHor == 1 && lenYVer/lenCbCrVer == 2
    subimg = [4 2 2];
elseif lenYHor/lenCbCrHor == 2 && lenYVer/lenCbCrVer == 2
    subimg = [4 2 0];
else
    error("Invalid Subsampling Vector.")
end

imgRec = convert2rgb(YCbCr{1},YCbCr{3},YCbCr{2}, subimg);

end