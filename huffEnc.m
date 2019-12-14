% Konstantinos Letros 8851
% Multimedia Systems Project
% Huffman Encoding

function huffStream = huffEnc(runSymbols)

global Nx
binaryStream = [];

for k = 1:size(runSymbols,1)
    
    % Find the Category of the non-zero DCT coeffs
    category = 0;
    
    if runSymbols(k,2)~=0
        while category<11
            if 2^(category-1)<=abs(runSymbols(k,2)) && abs(runSymbols(k,2))<=(2^category)-1
                break
            end
            category = category + 1;
        end
    else
        category = 0;
    end
    
    % Convert non-zero DCT Coeffs to binary
    
    if runSymbols(k,2)<0 % 1's Complement
        magnitude = double(~de2bi(-runSymbols(k,2),'left-msb'));
    else
        magnitude = de2bi(runSymbols(k,2),'left-msb');
    end
    
    
    if k==1
        % Find Index of Category in DC table
        huffmanInd = category + 1;
        DCvector=[Nx.DCTable{huffmanInd}-'0',magnitude];
        binaryStream = [binaryStream,DCvector];
    else
        % Find Index of Run/Category in AC table
        huffmanInd = 10*runSymbols(k,1)+(category+1);
        if runSymbols(k,1)==15 % index 152 --> ZRL
            huffmanInd = huffmanInd + 1;
        elseif isequal(runSymbols(k,:),[0 0]) % EOB
            magnitude = [];
        end
        ACvector = [Nx.ACTable{huffmanInd}-'0',magnitude];
        binaryStream = [binaryStream,ACvector];
    end
    
end

% Zero padding to create bytes
while (mod(length(binaryStream),8)) ~= 0
    binaryStream(end+1)=0;
end

% Convert Stream to bytes
i = 0;
k = 1;
huffStream = zeros(1,length(binaryStream)/8);
while i < length(binaryStream)
    huffStream(k) = uint8(bi2de(binaryStream(i+1:i+8),'left-msb'));
    i = i + 8;
    k = k + 1;
end

end