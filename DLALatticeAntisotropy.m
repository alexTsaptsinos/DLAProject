%% INVESTIGATING ANTISOTROPY IN DLA

clear
clf
clc

radius = input('radius: ');
tic
numberOfReps = 300;

cumulativeDensities = zeros(4*radius);
cumulativeNeighbourDensities = zeros(4*radius);
fractalDimensions = zeros(1,numberOfReps);
numberOfParticles = zeros(1,numberOfReps);
diameters = zeros(1,numberOfReps);
radii = zeros(1,numberOfReps);
edges = 0:pi/180:2*pi;
cumulativeRoundedAngles = zeros(1,360);
cos4N = 0;

for i = 1:numberOfReps
    fprintf(2,['ITERATE ' num2str(i) '\n']);
    
    [matrix,neighbourCountMatrix,particleAngles,fractalDimension...
        ,particleNumber,diameter] = DLALatticeMeakin(radius);
    
    cumulativeDensities = cumulativeDensities + matrix;
    cumulativeNeighbourDensities = cumulativeNeighbourDensities + neighbourCountMatrix;
    fractalDimensions(1,i) = fractalDimension;
    numberOfParticles(1,i) = particleNumber;
    radii(1,i) = diameter/2;
    
    % calculate <cos4x> quantative measure of anisotropy (from Alves/Ferreira
    % 2006
    numberOfAngles = size(particleAngles);
    numberOfAngles = numberOfAngles(1);
    cosAngles = cos(4*particleAngles);
    sumAngles = sum(cosAngles)/numberOfAngles;
    cos4N = cos4N + sumAngles;
    
    roundedAngles = histcounts(particleAngles,edges);
    cumulativeRoundedAngles = cumulativeRoundedAngles + roundedAngles;
    
end

fprintf(2,'FINISHED \n');

meanRadius = mean(radii);
disp(['Average Radius: ' num2str(meanRadius)]);
meanNumberOfParticles = mean(numberOfParticles);
disp(['Average Number Of Particles: ' num2str(meanNumberOfParticles)]);
meanFractalDimension = mean(fractalDimensions);
disp(['Average Fractal Dimension: ' num2str(meanFractalDimension)]);
disp(['<cos4x> quantitive anisotropy measure: ' num2str(cos4N)]);

clf
% Plots whole of all clusters on top of each other
figure(1)
imagesc(cumulativeDensities)
colormap(jet)
xlim([radius/2,7*radius/2])
ylim([radius/2,7*radius/2])
title(['Spatial Distribution average over ' num2str(numberOfReps) ' clusters'])

% Plots whole of all clusters neighbour distributions on top of each other
figure(2)
imagesc(cumulativeNeighbourDensities)
colormap(jet)
xlim([radius/2,7*radius/2])
ylim([radius/2,7*radius/2])
title(['Spatial Neighbour Distributions average over ' num2str(numberOfReps) ' clusters'])

% Plots only top right hand corner of all clusters on top of each other
figure(3)
imagesc(cumulativeDensities)
colormap(jet)
xlim([2*radius,7*radius/2])
ylim([radius/2,2*radius])
title(['Spatial Distribution average over ' num2str(numberOfReps) ' clusters'])

% Plots only top right hand corner of all clusters neighbour distributions
% on top of each other
figure(4)
imagesc(cumulativeNeighbourDensities)
colormap(jet)
xlim([2*radius,7*radius/2])
ylim([radius/2,2*radius])
title(['Spatial Neighbour Distributions average over ' num2str(numberOfReps) ' clusters'])

% Now plot graph of distribution of particle angles with cubic of best fit
first90CumulativeAngles = cumulativeRoundedAngles(1:90);
p = polyfit(0.5:1:89.5,first90CumulativeAngles,3);
a = 0.5:1:89.5;
b = polyval(p,a);
figure(5)
plot(a,b)
hold on
plot(first90CumulativeAngles,'s')
hold off

timeElapsed = toc;
disp(['Time Elapsed: ' num2str(timeElapsed)]);
