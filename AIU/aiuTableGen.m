clear all
close all
clc


d = 2.075;
numEls = 8;
nbits = 10;
delayVector = [16 32 63 125 250 500 1000 2000 4000 8000]*1e-3; %ns

minTheta = -60/180*pi;
maxTheta = 60/180*pi;
thetaPoints = 61;
c =  299792458;
d = d*.0254;

%% SSA
% This is a classic integer linear programming problem
theta = linspace(0,maxTheta,thetaPoints);
intcon = 1:nbits; %integer constraints
spacings = [8 16]*d;


lb = zeros(nbits,1);
ub = ones(nbits,1);
margin = delayVector(1);
f = delayVector;
A = repmat(delayVector,3,1);    
A(2:3,:)= A(2:3,:)*-1;

for j = 1:length(spacings) 
    dt = spacings(j)*sin(theta)/c/1e-9;
    for i = 1:length(dt)
        b = [dt(i)+margin;-dt(i)+margin;0];
        [x(i,j,:), fval(i,j)] = intlinprog(f,intcon,A,b,[],[],lb,ub);      
    end
end

%%
% subplot(211)
stem(squeeze(fval(:,1)))
% subplot(212)
hold on
stem(squeeze(fval(:,2)))

%% FA
% This is a classic integer linear programming problem
theta = linspace(0,maxTheta,thetaPoints);
intcon = 1:nbits; %integer constraints
spacings = [1:11]*8*d;

lb = zeros(nbits,1);
ub = ones(nbits,1);
margin = delayVector(1);
f = delayVector;
A = repmat(delayVector,3,1);    
A(2:3,:)= A(2:3,:)*-1;

%% 
for j = 1:length(spacings) 
    dt = spacings(j)*sin(theta)/c/1e-9;
    for i = 1:length(dt)
        b = [dt(i)+margin;-dt(i)+margin;0];
        [x(i,j,:), fval(i,j)] = intlinprog(f,intcon,A,b,[],[],lb,ub);      
    end
end

%%
for i = 1:length(spacings)
    hold on
    stem(-60:1:60,[flipud(fval(:,i))' fval(2:end,i)']');
end


