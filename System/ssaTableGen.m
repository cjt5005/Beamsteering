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

numStates = 3*(36+10);
numHWSAs = 3
intcon = 1:numStates; %integer constraints

spacings = [numHWSAs*numEls-1:-1:1]*d;
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

eSSA = [16 32 63 125 250 500 1000 2000 4000 8000]

eBlock = [repmat(eSSA,numEls,1) eHWSA];

e = blkdiag(eBlock,eBlock,eBlock);
%	[[eSSA1 eHWSA] [zeros]       [zeros]      ]
%   [[zeros]       [eSSA2 eHWSA] [zeros]      ]
%   [[zeros]       [zeros]       [eSSA3 eHWSA]]



eRelativeReference = (e - repmat(e(end,:),size(e,1),1));
eRelativeReference = eRelativeReference(1:end-1,:);
A = [eRelativeReference;-eRelativeReference];



%%
h = waitbar(0,'Reticulating Splines');
options=optimoptions(@intlinprog,'MaxTime',5,'IPPreprocess','advanced','CutGeneration','advanced');
xFullSSA = zeros(thetaPoints,numStates);
fvalFullSSA = zeros(1,thetaPoints);

for i = 1:length(dt)    

    tempMargin = 0;
    exitflag = 0;
    
    while exitflag <=0
        b = [dt(:,i) + repmat(tempMargin,length(spacings),1);
        -dt(:,i) + repmat(tempMargin,length(spacings),1)];
        [tempx,tempf , exitflag] = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
        tempMargin = tempMargin+1;
    end
        xFullSSA(i,:) = tempx;
        fvalFullSSA(i) = tempf;
    waitbar(i/length(dt))
end

%%
electronicsRelDelayVsAngle = eRelativeReference*xFullSSA';
stem(degTheta,electronicsRelDelayVsAngle')
