% Konstantinos Letros 8851
% Multimedia Systems Project
% JPEG Decode

function imgRec = JPEGdecode(JPEGenc, subimg, qScale)

YCbCr = cell(1,3);
types = ["Y","Cb","Cr"];

for count = 2 : length(JPEGenc)
    
    blockType = JPEGenc{count}.blkType;
    indHor  = JPEGenc{count}.indHor;
    indVer  = JPEGenc{count}.indVer;
    
    type = find(blockType==types);
    
    if type == 1
        qTable = JPEGenc{1}.qTableL;
    else
        qTable = JPEGenc{1}.qTableC;
    end
    
    if indHor == 1 && indVer == 1
        DCpred = 0;
    end
    
    runSymbols = huffDec(JPEGenc{count}.huffStream,type);
    
    qBlock = irunLength(runSymbols,DCpred);
    
    DCpred = qBlock(1,1);
    
    dctBlock = dequantizeJPEG(qBlock,qTable,qScale);
    block = iBlockDCT(dctBlock);
    
    YCbCr{type}((indHor-1)*8 + 1:indHor*8,(indVer-1)*8 + 1:indVer*8) = block;
    
end

imgRec = convert2rgb(YCbCr{1},YCbCr{2},YCbCr{3}, subimg);

end