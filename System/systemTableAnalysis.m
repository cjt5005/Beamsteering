close all
clear all
clc

load 'systemTableVariablesRev1'

%% Initialize variables
c =  299792458;
d = 2.05;
d = d*.0254;

theta = degTheta/180*pi;
numEls = 8;
%% HWSA Error Analysis

% Spacings between elements
spacings = [numEls-1:-1:1]*d;
% Time differences (pairs of element 1-7 and element 8) needed to steer beam at angle theta
dt = spacings'*sin(theta)/c/1e-12;

% Time delay differences (for each state variable) between each element and element 8
eRelativeHWSA = eHWSA - repmat(eHWSA(end,:),size(eHWSA,1),1);
eRelativeHWSA = eRelativeHWSA(1:end-1,:);

% Time differences(pairs of element 1-7 and element 8) derived from
% beamsteering algorithm
tHWSA = eHWSA*xElectronics';
dtActual = eRelativeHWSA*xElectronics';

stem(degTheta,dtActual')
hold on
plot(degTheta,dt')

figure()
for i = 1:numEls-1
    subplot(3,3,i)
    stem(degTheta,dtActual(i,:))
    hold on
    plot(degTheta,dt(i,:))
end

figure()
for i = 1:numEls-1
    subplot(3,3,i)
    stem(degTheta,dt(i,:)-dtActual(i,:))
end

err = (dt - dtActual)';

m = mean(err)
s = sqrt(var(err))
%% SSA Error Analysis

numHWSAs = 3;

% Spacings between HWSAs
spacings = (numEls*numHWSAs-1:-1:1)*d;
% Time differences (pairs of HWSAs 1-2 and HWSA 3) needed to steer beam at angle theta
dt = spacings'*sin(theta)/c/1e-12;

% Time delay for each HWSA
tSSA = eSSA*xAiuSSA';

for i = 1:numHWSAs
   tActual((i-1)*numEls+1:i*numEls,:) = repmat(tSSA(i,:),numEls,1) + tHWSA;    
end

% Time delay differences (for each state variable) between each element and element 8
dtActual = tActual - repmat(tActual(end,:),size(tActual,1),1);
dtActual = dtActual(1:end-1,:);

plot(degTheta,dt')
hold on
stem(degTheta,dtActual')

err = (dt - dtActual)';

figure()
stem(err)
m = mean(err)
s = sqrt(var(err))


%% SSA Error Analysis

numHWSAs = 12;

% Spacings between HWSAs
spacings = (numEls*numHWSAs-1:-1:1)*d;
% Time differences (pairs of HWSAs 1-2 and HWSA 3) needed to steer beam at angle theta
dt = spacings'*sin(theta)/c/1e-12;

% Time delay for each HWSA
tFA = eFA*xAiuFA';

for i = 1:numHWSAs
   tActual((i-1)*numEls+1:i*numEls,:) = repmat(tFA(i,:),numEls,1) + tHWSA;    
end

% Time delay differences (for each state variable) between each element and element 8
dtActual = tActual - repmat(tActual(end,:),size(tActual,1),1);
dtActual = dtActual(1:end-1,:);

plot(degTheta,dt')
hold on
stem(degTheta,dtActual')

err = (dt - dtActual)';

figure()
stem(err)
m = mean(err)
s = sqrt(var(err))
