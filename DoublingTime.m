%% Calculate Doubling Time

clear
clc
clf

% set up variables

radius = input('radius to double: ');

tic
particles = DLARectangular(radius);
% A_n has this many particles with r_n = radius
radiusTime = toc;
doubleParticles = DLARectangular(2*radius);
doubleRadiusTime = toc;
T = doubleParticles - particles;

beta = log(T)/log(particles)
