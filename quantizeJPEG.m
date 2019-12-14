% Konstantinos Letros 8851
% Multimedia Systems Project
% Quantizer

function qBlock = quantizeJPEG(dctBlock, qTable, qScale)

qBlock = round(dctBlock./(qScale*qTable));

end