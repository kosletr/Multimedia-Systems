%% Quantizer
function qBlock = quantizeJPEG(dctBlock, qTable, qScale)

qBlock = round(dctBlock./(qScale*qTable));

end