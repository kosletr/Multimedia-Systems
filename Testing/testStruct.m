load Testing/KyriaFile.mat
JPEGencKyr = JPEGenc;
load Testing/img.mat
JPEGencLet = JPEGenc;
clear JPEGenc
Tables;

score = 0;
i=2;
% for i = 2:length(JPEGencKyr)
    score = score  + isequal(JPEGencKyr{i}.huffStream,uint8(JPEGencLet{i}.huffStream));
% end
fprintf('Score: %f%% \n', 100*score/length(JPEGencLet))

runSymbolsKyr = huffDec(double(JPEGencKyr{i}.huffStream));
runSymbolsLet = huffDec(double(JPEGencLet{i}.huffStream));

qBlockKyr = irunLength(runSymbolsKyr,0);
qBlockLet = irunLength(runSymbolsLet,0);

dctBlockKyr = dequantizeJPEG(qBlockKyr,qTable{1},0.6);
dctBlockLet = dequantizeJPEG(qBlockLet,qTable{1},0.6);

blockKyr = uint8(iBlockDCT(dctBlockKyr))
blockLet = uint8(iBlockDCT(dctBlockLet))
