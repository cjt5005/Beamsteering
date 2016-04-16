clear all
close all
clc


d = 2.05;
numHWSAs = 12;
numEls = 8;
numStates = 46*numHWSAs;
electronicsDelayVector = [repmat([120 60 30],1,8) repmat([282 141],1,4) 542 272 542 272];
aiuDelayVector = [16 32 63 125 250 500 1000 2000 4000 8000]; %ps

delayVector = repmat([electronicsDelayVector aiuDelayVector],1,numHWSAs);

minTheta = -60/180*pi;
maxTheta = 60/180*pi;
thetaPoints = 61;
theta = linspace(minTheta,maxTheta,thetaPoints);
intcon = 1:numStates; %integer constraints

c =  299792458;
d = d*.0254;

%% Array Stuffs
% This is a classic integer linear programming problem

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
% row 1 => 30ps
% row 2 => 60ps
% row 3 => 120ps
% row 4 => 142ps
% row 5 => 282 ps
% row 6 => 272ps
% row 7 => 542ps
%% Full Array Diagram
%    13  	12  	11  	 10  	 9  	 8        7  	  6  	  5       4  	  3       2 
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
spacings = (numHWSAs*numEls-1:-1:1)*d;
startMargin = 15;
lb = zeros(numStates,1);
ub = ones(numStates,1);
f = ones(numStates,1);

% A Matrix
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
dt = spacings'*sin(theta)/c/1e-12;

%%
options=optimoptions(@intlinprog,'MaxTime',5,'IPPreprocess','advanced','CutGeneration','advanced','Display','final');
x = zeros(thetaPoints,numStates);
fval = zeros(1,thetaPoints);
parpool(4)
parfor i = 1:size(dt,2)
    
    percent_complete = i/size(dt,2)*100
    tempMargin = startMargin;
    exitflag = 0;
    
    while exitflag <=0
        b = [dt(:,i) + repmat(tempMargin,length(spacings),1);
        -dt(:,i) + repmat(tempMargin,length(spacings),1)];
        [tempx,tempf , exitflag] = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
        tempMargin = tempMargin+1;
    end
        x(i,:) = tempx;
        fval(i) = tempf;
end

%%

for i = 1:numHWSAs*numEls-1
    stem(linspace(minTheta/pi*180,maxTheta/pi*180,thetaPoints),(e(i,:)-e(end,:))*x')
    hold on
end

delete(gcp('nocreate'))
    