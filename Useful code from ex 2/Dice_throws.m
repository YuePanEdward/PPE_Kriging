% Part 2.2 d covariance estimation

%% Clean up
clear all;
close all;
clc;

%% Dice throws
D = [5, 4, 2, 4, 5, 3, 4, 1, 1, 4, 2, 5, 4, 1, 2, 5, 1, 2, 2, 6, 6, 3, 6, 4, 6, 4, 1, 3, 3, 2];

%% Calculating Xi
X = zeros(1,length(D));

for i = 1: length(X)
    if i == 1
        X(i) = 1.75 + 0.5*D(i);
        
    else
        X(i) = 0.5*D(i-1) + 0.5*D(i);
        
    end
    
end

%% Posible distances 
p = [0:length(X)-1];

%% Calculating the semi-variogram in relation to the distance
semi = zeros(0,length(p));

for i = 1:length(p)
    semi(i) = semiVariogram(p(i),X);
end

%% Calculating the experimental variance
meanX = mean(X);
difVar = X-mean(X);
squareVar = difVar .* difVar;
sumSquareVar = sum(squareVar);
Var = 1/length(X)*sumSquareVar;

%% Calculating cov(x_t,x_t+p)
cov = Var - semi;

%% Plotting the cov
figure;
pFigure = [-flip(p(2:end)), p];
covFigure = [flip(cov(2:end)) cov];
plot(pFigure,covFigure,'*-');
title('Coputed covariances in relationship to p');
xlabel('p-Values');
ylabel('Value of the covariance');

%% Export
f = gcf;
exportgraphics(f,'covPlotDices.png')
