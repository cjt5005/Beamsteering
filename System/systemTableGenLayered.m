clear all
close all
clc

%% Antenna Electronics Diagram
%    1  	   2  	   3	  4  	     5  	    6  	    7       8
%   \|/ 	  \|/     \|/ 	 \|/ 	    \|/ 	   \|/     \|/ 	   \|/
%    |  	   |  	|  	  |  	     |  	    |       |  	    |
%   TD1     TD4     TD7     TD10     TD13     TD16     TD19     TD22
%    |       |       |       |         |        |        |        |
%   TD2     TD5     TD8     TD11     TD14     TD17     TD20     TD23
%    |       |       |       |         |        |        |        |
%   TD3     TD6     TD9     TD12     TD15     TD18     TD21     TD24
%	 |_______|       |_______|         |________|        |________|
%	    |                |          	      |           	  |
%	   TD25            TD27                TD29              TD31
%	    |           	    |          	      |           	  |
%	   TD26            TD28                TD30              TD32
%        |_______________|                   |_________________|
%                   |                                      |
%			   TD33                 		       TD35
%			    |                  			   |
%			   TD34                 		       TD36
%			    |______________________________________|
%								 |
%		 					  RF Input
%
% row 1 => 30ps
% row 2 => 60ps
% row 3 => 120ps
% row 4 => 142ps
% row 5 => 282 ps
% row 6 => 272ps
% row 7 => 542ps
%% SSA Diagram
%    1  	 2  	 3  	 
%   \|/     \|/     \|/ 	 
%    |  	 |  	 |  	
%   TD1    TD11    TD21    
%    |       |       |     
%   TD2    TD12    TD22     
%    |       |       |     
%   ...     ...     ...
%    |       |       |
%   TD10    TD20   TD30     
%	 |_______|_______|
%            |               
%	         |            
%	         |           	
%	     RF Input          
%            
%
% row 1 => 16ps
% row 2 => 32ps
% row 3 => 63ps
% row 4 => 125ps
% row 5 => 250ps
% row 6 => 500ps
% row 7 => 1000ps
% row 8 => 2000ps
% row 9 => 4000ps
% row 10 => 8000ps
%% Full Array Diagram
%    1  	 2  	 3  	 4  	 5  	 6        7  	  8  	  9       10  	  11	  12  
%   \|/     \|/     \|/ 	\|/     \|/     \|/      \|/     \|/     \|/      \|/     \|/     \|/ 
%    |  	 |  	 |  	 |  	 |  	 |        |  	  |  	  |        |  	   |  	   |  
%   TD1    TD11    TD21     TD31    TD41    TD51    TD61    TD71     TD81    TD91    TD101   TD111  
%    |       |       |       |       |       |        |       |       |        |       |       |   
%   TD2    TD12    TD22     TD32    TD42    TD52    TD62    TD72     TD82    TD92    TD102   TD112  
%    |       |       |       |       |       |        |       |       |        |       |       |   
%   ...     ...     ...     ...     ...     ...      ...     ...     ...      ...     ...     ...
%    |       |       |       |       |       |        |       |       |        |       |       |
%   TD10    TD20   TD30     TD40    TD50   TD60     TD70    TD80    TD90     TD100   TD110   TD120  
%	 |_______|_______|_______|_______|_______|________|_______|_______|________|_______|_______|
%            									 |                                                       
%	         									 |                                                       
%	         									 |           	                                          
%	     									  RF Input              
% row 1 => 16ps
% row 2 => 32ps
% row 3 => 63ps
% row 4 => 125ps
% row 5 => 250ps
% row 6 => 500ps
% row 7 => 1000ps
% row 8 => 2000ps
% row 9 => 4000ps
% row 10 => 8000ps

%%
c =  299792458;
d = 2.05;
d = d*.0254;

numEls = 8;
electronicsDelayVector = [repmat([120 60 30],1,8) repmat([282 141],1,4) 542 272 542 272];
aiuDelayVector = [16 32 63 125 250 500 1000 2000 4000 8000]; %ps

minTheta = -60/180*pi;
maxTheta = 60/180*pi;
thetaPoints = 121;
theta = linspace(minTheta,maxTheta,thetaPoints);
degTheta = theta*180/pi;


%% HWSA
% This is a classic integer linear programming problem

numStates = 36;
intcon = 1:numStates; %integer constraints

spacings = [numEls-1:-1:1]*d;
dt = spacings'*sin(theta)/c/1e-12;

lb = zeros(numStates,1);
ub = ones(numStates,1);
f = ones(numStates,1);

% A Matrix
eHWSA =   [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0;...
           0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0;...
           0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1;...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 0 0 1 1].*repmat(electronicsDelayVector,numEls,1);

eRelativeReference = (eHWSA - repmat(eHWSA(end,:),size(eHWSA,1),1));
eRelativeReference = eRelativeReference(1:end-1,:);
A = [eRelativeReference;-eRelativeReference];



%%
options=optimoptions(@intlinprog,'display','off');
xElectronics = zeros(thetaPoints,numStates);
fvalElectronics = zeros(1,thetaPoints);

for i = 1:length(dt)    
    percent_complete = i/length(dt)*100
    tempMargin = 0;
    exitflag = 0;
    
    while exitflag <=0
        b = [dt(:,i) + repmat(tempMargin,length(spacings),1);
        -dt(:,i) + repmat(tempMargin,length(spacings),1)];
        [tempx,tempf , exitflag] = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
        tempMargin = tempMargin+1;
    end
        xElectronics(i,:) = tempx;
        fvalElectronics(i) = tempf;
end

%%
electronicsRelDelayVsAngle = eRelativeReference*xElectronics';
stem(degTheta,electronicsRelDelayVsAngle')

%% SSA
% This is a classic integer linear programming problem

numStates = 30;
numEls = 8;
numHWSAs = 3;
intcon = 1:numStates; %integer constraints

spacings = [numHWSAs-1:-1:1]*d*numEls;
dt = spacings'*sin(theta)/c/1e-12;

lb = zeros(numStates,1);
ub = ones(numStates,1);
f = ones(numStates,1);

% A Matrix
eSSA =   blkdiag(aiuDelayVector,aiuDelayVector,aiuDelayVector);

eRelativeReference = (eSSA - repmat(eSSA(end,:),size(eSSA,1),1));
eRelativeReference = eRelativeReference(1:end-1,:);
A = [eRelativeReference;-eRelativeReference];

%%
options=optimoptions(@intlinprog,'display','off');
xAiuSSA = zeros(thetaPoints,numStates);
fAiuSSA= zeros(1,thetaPoints);

for i = 1:length(dt)    
    percent_complete = i/length(dt)*100
    tempMargin = 0;
    exitflag = 0;
    
    while exitflag <=0
        b = [dt(:,i) + repmat(tempMargin,length(spacings),1);
        -dt(:,i) + repmat(tempMargin,length(spacings),1)];
        [tempx,tempf , exitflag] = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
        tempMargin = tempMargin+1;
    end
        xAiuSSA(i,:) = tempx;
        fAiuSSA(i) = tempf;
end

%%
aiuSsaRelDelayVsAngle = eRelativeReference*xAiuSSA';
stem(degTheta,aiuSsaRelDelayVsAngle')

%% FA
% This is a classic integer linear programming problem

numStates = 120;
numEls = 8;
numHWSAs = 12;
intcon = 1:numStates; %integer constraints

spacings = [numHWSAs-1:-1:1]*d*numEls;
dt = spacings'*sin(theta)/c/1e-12;

lb = zeros(numStates,1);
ub = ones(numStates,1);
f = ones(numStates,1);

% A Matrix
eFA =   blkdiag(aiuDelayVector,aiuDelayVector,aiuDelayVector,...
                aiuDelayVector,aiuDelayVector,aiuDelayVector,...
                aiuDelayVector,aiuDelayVector,aiuDelayVector,...
                aiuDelayVector,aiuDelayVector,aiuDelayVector);

eRelativeReference = (eFA - repmat(eFA(end,:),size(eFA,1),1));
eRelativeReference = eRelativeReference(1:end-1,:);
A = [eRelativeReference;-eRelativeReference];

%%
options=optimoptions(@intlinprog,'display','off');
xAiuFA = zeros(thetaPoints,numStates);
fAiuFA= zeros(1,thetaPoints);

for i = 1:length(dt)    
    percent_complete = i/length(dt)*100
    tempMargin = 0;
    exitflag = 0;
    
    while exitflag <=0
        b = [dt(:,i) + repmat(tempMargin,length(spacings),1);
        -dt(:,i) + repmat(tempMargin,length(spacings),1)];
        [tempx,tempf , exitflag] = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
        tempMargin = tempMargin+1;
    end
        xAiuFA(i,:) = tempx;
        fAiuFA(i) = tempf;
end

%%
aiuFaRelDelayVsAngle = eRelativeReference*xAiuFA';
stem(degTheta,aiuFaRelDelayVsAngle')