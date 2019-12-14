% Konstantinos Letros 8851
% Multimedia Systems Project
% De-Quantizer

function dctBlock = dequantizeJPEG(qBlock, qTable, qScale)

dctBlock = qBlock.*(qScale*qTable);

end