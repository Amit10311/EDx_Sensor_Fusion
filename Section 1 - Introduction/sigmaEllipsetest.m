%% %Here we give parameters for a Gaussian density. The parameter mu is the mean, and P is the covariance matrix.

mu = [-2; 1];
P = [4, -2; -2 2];

%Call your function.
xy = sigmaEllipse2D(mu, P);

%Now plot the generated points. 
%You should see an elongated ellipse stretching from the top left corner to the bottom right. 
figure(1);
h1 = plot(xy(1,:), xy(2,:));

%Set the scale of x and y axis to be the same. 
% This should be done if the two variables are in the same domain,e.g. both are measured in meters.
axis equal
hold on
title("Elongated ellipse ")
%Also plot a star where the mean is, and make it have the same color as the ellipse.
plot(mu(1), mu(2), '*', 'color', h1.Color);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%Pretest

%There should be two rows in the output (Pretest)
% randomize a 2D Gaussian distribution

mu = 5*randn(2,1);
s = 2*rand(2);
P = s'*s;
xy = sigmaEllipse2D(mu, P)

assert(size(xy,1) == 2, 'There should be two rows in the output');

 %   2. Assuming roughly equidistant points on the ellipse, is the mean of the ellipse correct? 
% randomize a 2D Gaussian distribution
mu = 5*randn(2,1);
s = 2*rand(2);
P = s'*s;

xy = sigmaEllipse2D(mu, P)

xy = unique(xy', 'rows')'

assert(all(abs(mean(xy, 2) - mu)./norm(mu) < 1e-1), ...
    'The mean of the ellipse points should roughly coincide with the distribution mean');

%% %  3. Are the generated points on an ellipse of the right shape 

% randomize a 2D Gaussian distribution
mu = 5*randn(2,1);
s = 2*rand(2);
P = s'*s;
level = 1+rand*2;
npoints = 100+round(25*rand);

xy = sigmaEllipse2D(mu, P, level, npoints)

n = size(xy,2);
r2 = zeros(n,1);

for k = 1:n
    r2(k) = (xy(:,k)-mu)'/P*(xy(:,k)-mu);
end

mr2 = mean(r2);

assert(all(abs(r2-mr2) < 1e-5), 'The ellipse has the wrong shape') 

figure(1);
h1 = plot(xy(1,:), xy(2,:));
title("Generated points on an ellipse of the right shape ")
%Set the scale of x and y axis to be the same. 
% This should be done if the two variables are in the same domain,e.g. both are measured in meters.
axis equal
hold on

%Also plot a star where the mean is, and make it have the same color as the ellipse.
plot(mu(1), mu(2), '*', 'color', h1.Color);

hold off

%% % 4.  Check that the level option is used to control the size of the ellipse.

% randomize a 2D Gaussian distribution
mu = 5*randn(2,1);
s = 2*rand(2);
P = s'*s;
level = 1+rand*2;

xy = sigmaEllipse2D(mu, P, level)

n = size(xy,2);
r2 = zeros(n,1);
for k = 1:n
    r2(k) = (xy(:,k)-mu)'/P*(xy(:,k)-mu);
end

mr = sqrt(mean(r2))

assert(abs(mr-level) < 1e-3, ['The ellipse has the wrong size.' ...
    'The third parameter should control the size of the ellipse.'])

figure(1);
h1 = plot(xy(1,:), xy(2,:));

%Set the scale of x and y axis to be the same. 
% This should be done if the two variables are in the same domain,e.g. both are measured in meters.
axis equal
hold on

%Also plot a star where the mean is, and make it have the same color as the ellipse.
plot(mu(1), mu(2), '*', 'color', h1.Color);
hold off

%% 

function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )
%SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
% a sigma level in the Gaussian density defined by mean and covariance.
%
%Input:
%   MU          [2 x 1] Mean of the Gaussian density
%   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
%   LEVEL       Which sigma level curve to plot. Can take any positive value, 
%               but common choices are 1, 2 or 3. Default = 3.
%   NPOINTS     Number of points on the ellipse to generate. Default = 32.
%
%Output:
%   XY          [2 x npoints] matrix. First row holds x-coordinates, second
%               row holds the y-coordinates. First and last columns should 
%               be the same point, to create a closed curve.


%Setting default values, in case only mu and Sigma are specified.
if nargin < 3
    level = 3;
end
if nargin < 4
    npoints = 32;
end

%Your code here

% linspace of theta range
angle_samples = linspace(0.0, 2* pi, npoints);

polar_tf = [cos(angle_samples); sin(angle_samples)];

%using equation (2)
xy = mu + level * (sqrtm(Sigma) * polar_tf);
end



