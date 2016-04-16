function [ aiuStates,diuStates ] = exportStatesSSA(xSSA)

    block = [10 36];
    blockSize = sum(block);
    numBlocks = size(xSSA,2)/blockSize;
    
    for i = 1:numBlocks   
        aiuStates(:,:,i) = xSSA(:,(i-1)*blockSize+1:(i-1)*blockSize+block(1));
        diuStates(:,:,i) = xSSA(:,(i-1)*blockSize+block(1)+1:i*blockSize);    
    end

end

