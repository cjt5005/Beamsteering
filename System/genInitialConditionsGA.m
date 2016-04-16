    clear all
    close all
    clc
    
    load exampleSolution10degrees
    c =  299792458;
    theta = 10/180*pi;
    d = 2.05;
    numHWSAs = 12;
    numEls = 8;
    d = d*.0254;
    numStates = 46*numHWSAs;
    
    intcon = 1:numStates;
    ub = ones(1,numStates);
    lb = zeros(1,numStates);
    spacings = (numHWSAs*numEls-1:-1:1)*d;
    dt = spacings'*sin(theta)/c/1e-12;
            electronicsDelayVector = [repmat([120 60 30],1,8) repmat([282 141],1,4) 542 272 542 272];
        aiuDelayVector = [16 32 63 125 250 500 1000 2000 4000 8000];
    diuDelayMatrix = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0;...
           0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 0 0 1 1].*repmat(electronicsDelayVector,numEls,1);

aiuDelayMatrix = repmat(aiuDelayVector,numEls,1);
blockMat = [aiuDelayMatrix diuDelayMatrix];

e = blkdiag(blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat,blockMat);
% e = blkdiag(blockMat,blockMat,blockMat);

for i = 1:numHWSAs*numEls-1
    A(i,:) = e(i,:)-e(end,:);
end
 
% A = [A;-A;A2]; 
A = [A;-A];
    startMargin = 15;
    b = [dt + repmat(startMargin,length(spacings),1);
        -dt + repmat(startMargin,length(spacings),1)];

initialGuess = derp(:);

clearvars -except 'A' 'b' 'ub' 'lb' 'intcon' 'numStates' 'initialGuess'
