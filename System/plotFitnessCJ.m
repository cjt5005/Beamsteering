function state = plotFitnessCJ(~,state,flag)
    
    if  size(state.Score,2) > 1
        title('Best Individual Plot: not available','interp','none');
        return;
    end
    
    switch flag
        case 'init'
        return;
        
        case 'iter'


        x = state.Population(state.Score == min(state.Score),:)';
        if size(x,2) >1
            x = x(:,1);
        end
        
        c =  299792458;
        theta = 10/180*pi;
        d = 2.05;
        numHWSAs = 12;
        numEls = 8;
        d = d*.0254;



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


        spacings = (numHWSAs*numEls-1:-1:1)*d;
        dt = spacings'*sin(theta)/c/1e-12;
        

        calcT = e*x;
        calcT = calcT - calcT(end);
        calcT = calcT(1:end-1);

        

        plot(dt)
        hold on
        plot(calcT)
        
        hold off
    end

end

