function [matrix,neighbourCountMatrix,particleAngles,fractalDimension,medianNumberOfSteps,particleNumber,diameter] = DLALatticeBasic(radius)
%% DLA SIMULATION USING SQUARE MATRIX FOR ON LATTICE - ESCAPE RADIUS ONLY 2*INITIAL START

% Set up variables - what radius size do you want the cluster
% radius = input('radius: ');
tic

%% Initial Setup:
% We create a matrix of zeros and then plant a seed at the centre of the
% matrix. Change entry of matrix to 1 if it is aggregated.
%
% This is the matrix set up, we start from inner square which grows as the cluster grows and consider
% particle escaped if it reaches the outer square
%  -----------
% |    ----    |
% |   |    |   |
% |   |    |   |
% |    ----    |
%  -----------

matrix = zeros(4*radius);

middle = 2*radius;

matrix(middle,middle) = 1;

% Initialise our variables

particleNumber = 1; %integer
maximumDistance = 0;
endScript = 0; %bool
numberOfStepsMatrix = []; %matrix to record number of steps taken until aggregated

% We record the stuck particles
% in a matrix (stuck_particles) which records the x position, y position,
% distance from origin, angle from origin
stuck_particles = zeros(1,4);
stuck_particles(1,1) = middle;
stuck_particles(1,2) = middle;

%% Random Walk Script

while endScript == 0
    
    %first decide where the particle starts from
    R = maximumDistance + 5;
    randAngle = 2*pi*rand;
    
    if randAngle < pi/2
        x = R*sin(randAngle);
        x = round(x);
        y = R*cos(randAngle);
        y = round(y);
        x = middle + x;
        y = middle - y;
    elseif randAngle < pi
        theta = randAngle - pi/2;
        x = R*cos(theta);
        x = round(x);
        y = R*sin(theta);
        y = round(y);
        x = middle + x;
        y = middle + y;
    elseif randAngle < 3*pi/2
        theta = randAngle - pi;
        x = R*sin(theta);
        x = round(x);
        y = R*cos(theta);
        y = round(y);
        x = middle - x;
        y = middle + y;
    else
        theta = randAngle - 3*pi/2;
        x = R*cos(theta);
        x = round(x);
        y = R*sin(theta);
        y = round(y);
        x = middle - x;
        y = middle - y;
    end
    
    % Now we have a starting point on the edge of the inner square of the matrix. Now we want
    % to perform random walk in the matrix until we are next to an entry
    % which is 1. Set aggregate to false.
    
    aggregate = 0; %bool to exit below while loop
    escape = 0; %bool to exit below while loop
    numberOfSteps = 0; %record number of steps taken until aggregated
    
    while (aggregate == 0) && (escape == 0)
        numberOfSteps = numberOfSteps + 1;
        randWalk = rand;
        if randWalk < 1/4
            % go right
            x = x+1;
        elseif randWalk < 2/4
            % go left
            x = x-1;
        elseif randWalk < 3/4
            % go up
            y = y-1;
        else
            % go down
            y = y+1;
        end
        
        xdistanceFromCenter = abs(x - middle);
        ydistanceFromCenter = abs(y - middle);
        distanceFromCenter = sqrt(xdistanceFromCenter^2 + ydistanceFromCenter^2);
        
        % escape test
        %         if distanceFromCenter > 2*R
        %             escape = 1;
        %         else
        if (x == 1) || (x == 4*radius) || (y == 1) || (y == 4*radius)
            escape = 1;
        else
            % aggregate test, only if not escaped
            if (matrix(y,x+1) + matrix(y,x-1) + matrix(y+1,x) + matrix(y-1,x)) ~= 0
                aggregate = 1;
            end
        end
        
        %We keep repeating this until the particle is finally stuck or
        %escapes
        
    end
    
    %if escape
    %disp('particle escaped')
    %end
    
    if aggregate
        particleNumber = particleNumber + 1;
        stuck_particles(particleNumber,1) = x;
        stuck_particles(particleNumber,2) = y;
        stuck_particles(particleNumber,3) = distanceFromCenter;
        if distanceFromCenter > maximumDistance
            maximumDistance = distanceFromCenter;
        end
        % add angle to stuck_particles  matrix
        xquadrant = sign(x - middle);
        yquadrant = sign(y - middle);
        if xquadrant == 1
            if yquadrant == 1
                particleAngle = atan(ydistanceFromCenter/xdistanceFromCenter) + pi/2;
            elseif yquadrant ==  -1
                particleAngle = atan(xdistanceFromCenter/ydistanceFromCenter);
            else
                particleAngle = pi/2;
            end
        elseif xquadrant == -1
            if yquadrant == 1
                particleAngle = atan(xdistanceFromCenter/ydistanceFromCenter) + pi;
            elseif yquadrant == -1
                particleAngle = atan(ydistanceFromCenter/xdistanceFromCenter) + 3*pi/2;
            else
                particleAngle = 3*pi/2;
            end
        else
            if yquadrant == 1
                particleAngle = pi;
            else
                particleAngle = 0;
            end
        end
        stuck_particles(particleNumber,4) = particleAngle;
        
        disp(['particle aggregated: ' num2str(particleNumber)])
        numberOfStepsMatrix = [numberOfSteps; numberOfStepsMatrix];
        
        matrix(y,x)=particleNumber;
        
        % end of script check
        
        if distanceFromCenter > radius
            endScript = 1;
            % calculate fractal dimension
            
            diameter = 2*distanceFromCenter;
            fractalDimension = log(particleNumber)/log(distanceFromCenter);
            
        end
    end
    
end

%% Plot graph

matrix = matrix/particleNumber;
imagesc(matrix)
colormap(jet)
title(['DLA with ' num2str(particleNumber) ' particles'])
%text(width/3,3*R/2 + R/6,['Fractal Dimension: ' num2str(fractalDimension)]);
%text(R/3,3*R/2,['Radius: ' num2str(maximumDistance)]);
timeElapsed = toc;
%text(R/3,3*R/2 - R/6, ['Time Elapsed: ' num2str(timeElapsed)]);
%xlabel(num2str(2*width))
%ylabel(num2str(2*width))
%axis equal
xlim([radius,3*radius])
ylim([radius,3*radius])
axis off

%% Display Outputs

disp(['Number of particles: ' num2str(particleNumber)]);
disp(['Fractal Dimension: ' num2str(fractalDimension)]);
disp(['Diameter: ' num2str(diameter)]);
disp(['Time Elapsed: ' num2str(timeElapsed)]);
meanNumberOfSteps = mean(numberOfStepsMatrix);
medianNumberOfSteps = median(numberOfStepsMatrix);
%sdNumberOfSteps = std(numberOfStepsMatrix);
particleAngles = stuck_particles(:,4);
disp(['Mean Number Of Steps until aggregated: ' num2str(meanNumberOfSteps)]);
disp(['Median Number Of Steps until aggregated: ' num2str(medianNumberOfSteps)]);
%disp(['Standard Deviation of Number Of Steps until aggregated: ' num2str(sdNumberOfSteps)]);

neighbourCountMatrix = zeros(4*radius);
for i = 1:4*radius
    for j = 1:4*radius
        neighbourCount = 0;
        if j ~= 1
            if matrix(j-1,i) ~= 0
                neighbourCount = neighbourCount + 1;
            end
        end
        if j ~= 4*radius
            if matrix (j+1,i) ~= 0
                neighbourCount = neighbourCount + 1;
            end
        end
        if i ~= 1
            if matrix(j,i-1) ~= 0
                neighbourCount = neighbourCount + 1;
            end
        end
        if i ~= 4*radius
            if matrix (j,i+1) ~= 0
                neighbourCount = neighbourCount + 1;
            end
        end
        neighbourCountMatrix(j,i) = neighbourCount;
    end
    
end