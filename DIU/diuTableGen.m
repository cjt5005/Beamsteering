clear all
close all
clc
%    1  	   2  	3  	  4  	     5  	    6  	 7        8
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
%%
d = 2.05;
numStates = 36;
delayVector = [repmat([120 60 30],1,8) repmat([282 141],1,4) 542 272 542 272]; %ps


maxTheta = 90/180*pi;
thetaPoints = 91;
c =  299792458;
d = d*.0254;

%% HWSA
% This is a classic integer linear programming problem
theta = linspace(0,maxTheta,thetaPoints);
intcon = 1:numStates; %integer constraints
spacings = [7:-1:1]*d;

lb = zeros(numStates,1);
ub = ones(numStates,1);
f = ones(numStates,1);

% A Matrix
e1  = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0].*delayVector;
e2  = [0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0].*delayVector;
e3  = [0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0].*delayVector;
e4 = [0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0].*delayVector;
e5 = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1].*delayVector;
e6 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1].*delayVector;
e7 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1].*delayVector;
e8 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 0 0 1 1].*delayVector;

A = [e1-e8;
     e2-e8;
     e3-e8;
     e4-e8;
     e5-e8;
     e6-e8;
     e7-e8];
 
A = [A;-A]; 

dt = spacings'*sin(theta)/c/1e-12;

%%
options=optimoptions(@intlinprog,'display','off');
x = zeros(thetaPoints,numStates);
fval = zeros(1,thetaPoints);

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
        x(i,:) = tempx;
        fval(i) = tempf;
end

%%
s7val = (e1-e8)*x';
s8val =  (e2-e8)*x';
s9val  =  (e3-e8)*x';
s10val  =  (e4-e8)*x';
s11val  =  (e5-e8)*x';
s12val  =  (e6-e8)*x';
s13val  =  (e7-e8)*x';
%%

stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s7val) s7val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s8val) s8val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s9val) s9val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s10val) s10val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s11val) s11val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s12val) s12val(2:end)]');
hold on
stem(-thetaPoints+1:1:thetaPoints-1,[-fliplr(s13val) s13val(2:end)]');

