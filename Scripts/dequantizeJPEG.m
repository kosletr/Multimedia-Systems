%% De-Quantizer
function dctBlock = dequantizeJPEG(qBlock, qTable, qScale)

dctBlock = qBlock.*(qScale*qTable);

end

