load Testing/KyriaFile.mat
JPEGencKyr = JPEGenc;
load Testing/img.mat
JPEGencLet = JPEGenc;
clear JPEGenc
Tables;

i=10;
JPEGencKyr{2}.huffStream(1:i)
uint8(JPEGencLet{2}.huffStream(1:i))

score = 0;
 for i = 2:length(JPEGencKyr{2}.huffStream)
    score = score  +  isequal(double(JPEGencKyr{2}.huffStream(i)),(JPEGencLet{2}.huffStream(i)));
 end
fprintf('Score: %f%% \n', 100*score/length(JPEGencLet{2}.huffStream))


runSymbolsKyr = huffDec(double(JPEGencKyr{2}.huffStream));
runSymbolsLet = huffDec(double(JPEGencLet{2}.huffStream));

qBlockKyr = irunLength(runSymbolsKyr,0);
qBlockLet = irunLength(runSymbolsLet,0);

dctBlockKyr = dequantizeJPEG(qBlockKyr,qTable{1},0.6);
dctBlockLet = dequantizeJPEG(qBlockLet,qTable{1},0.6);

blockKyr = uint8(iBlockDCT(dctBlockKyr));
blockLet = uint8(iBlockDCT(dctBlockLet));
