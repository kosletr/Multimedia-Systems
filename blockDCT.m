% Konstantinos Letros 8851
% Multimedia Systems Project
% Discrete Cosine Transform of a block

function dctBlock = blockDCT(block)

dctBlock =  dct2(double(block)-128);

end