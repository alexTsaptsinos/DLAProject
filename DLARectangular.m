function [n] = DLARectangular(radius);

%% DLA SIMULATION USING SQUARE MATRIX


% Set up variables
%width = input('width: ');
tic
% numberOfParticles = input('Number of Particles: ');

%% Initial Setup:
% We create a matrix of zeros and then plant a seed at the centre of the
% matrix. Change entry of matrix to 1 if it is aggregated.
%
% This is the matrix set up, we start from inner square and consider
% particle escaped if it reaches the outer square
%  -----------
% |    ----    |
% |   |    |   |
% |   |    |   |
% |    ----    |
%  -----------

matrix = zeros(5*radius,5*radius);

middlex = floor(size(matrix,1)/2);
middley = floor(size(matrix,1)/2);

matrix(middley,middlex) = 0.1;

% We want to let a particle randomly enter from the edges of the matrix so
% first let's choose a side

particleNumber = 1; %integer
x = 0; %integer
y = 0; %integer
aggregate = 0; %bool
endScript = 0; %bool
halfRadius = 0; %bool
escape = 0; %bool
startingRadius = 5; %double

%% Random Walk Script

while endScript == 0
    
    %first decide where the particle starts from    
    randAngle = 2*pi*rand;
    hold on
    
    if randAngle < pi/2
        x = round(middlex + startingRadius*sin(randAngle));
        y = round(middley + startingRadius*cos(randAngle));
    elseif randAngle < pi
        theta = randAngle - pi/2;
        x = round(middlex + startingRadius*cos(theta));
        y = round(middley - startingRadius*sin(theta));
    elseif randAngle < 3*pi/2
        theta = randAngle - pi;
        x = round(middlex - startingRadius*sin(theta));
        y = round(middley - startingRadius*cos(theta));
    else
        theta = randAngle - 3*pi/2;
        x = round(middlex - startingRadius*cos(theta));
        y = round(middley + startingRadius*sin(theta));
    end 
    
    % Now we have a starting point on the edge of the inner square of the matrix. Now we want
    % to perform random walk in the matrix until we are next to an entry
    % which is 1. Set aggregate to false.
    
    aggregate = 0;
    escape = 0;
    
    while (aggregate == 0) && (escape == 0)
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
        
        % escape test

        distanceFromCenter = sqrt((x-middlex)^2 + (y - middley)^2);
        if distanceFromCenter > 2*startingRadius
            escape = 1;
        end
        
        % aggregate test, only if not escaped
        
        if escape == 0
            if (matrix(y,x+1) + matrix(y,x-1) + matrix(y+1,x) + matrix(y-1,x)) ~= 0
                aggregate = 1;
            end
        end
        
        %We keep repeating this until the particle is finally stuck
        
    end
    
    if escape
        disp('particle escaped')
    end
    
    if aggregate
        particleNumber = particleNumber + 1;
        xdif = x - middlex;
        ydif = y - middley;
        distanceFromCenter = sqrt(ydif^2 + xdif^2);
        if distanceFromCenter > startingRadius
            startingRadius = distanceFromCenter + 5;
        end
        % half radius check
        if halfRadius == 0;
            if distanceFromCenter >= radius/2
                disp('HALF RADIUS');
                halfRadius = 1;
                halfParticleNumber = particleNumber;
            end
        end
        if distanceFromCenter >= radius
            endScript = 1;
            fractdim = log(particleNumber)/log(distanceFromCenter);
            doubleTime = log(particleNumber - halfParticleNumber)/log(particleNumber);
        end
        disp(['particle aggregated: ' num2str(particleNumber)])
        if particleNumber < 500
            matrix(y,x) = 0.1;
        elseif particleNumber < 1000
            matrix(y,x) = 0.2;
        elseif particleNumber < 1500
            matrix(y,x) = 0.3;
        elseif particleNumber < 2000
            matrix(y,x) = 0.4;
        elseif particleNumber < 2500
            matrix(y,x) = 0.5;
        elseif particleNumber < 3000
            matrix(y,x) = 0.6;
        elseif particleNumber < 3500
            matrix(y,x) = 0.7;
        elseif particleNumber < 4000
            matrix(y,x) = 0.8;
        elseif particleNumber < 4500
            matrix(y,x) = 0.9;
        else
            matrix(y,x) = 1.0;
        end
        
        % end of script check
        
%         if (x==ceil(width/2) || x==ceil(3*width/2) || y==ceil(3*width/2) || y==ceil(width/2))
%             endScript = 1;
%             % now want to calculate diameter of fractal
%             % we could change to make sure that we always add particles until we hit
%             % the edge, then diameter we could work out
%             fractdim = log(particleNumber)/log(radius);
%             
%         end
    end
    
end

%Tell why escaped, if true
% if (endScript == 1)
%     disp ('Script ended early as aggregation was about to exceed matrix dimensions');
%
% end

% uncomment this for diagram to just appear at end

%% Plot graph

figure(1)
imagesc(matrix)
colormap(jet)
title(['DLA with ' num2str(particleNumber) ' particles'])
%text(width/3,3*R/2 + R/6,['Fractal Dimension: ' num2str(fractalDimension)]);
%text(R/3,3*R/2,['Radius: ' num2str(maximumDistance)]);
timeElapsed = toc;
%text(R/3,3*R/2 - R/6, ['Time Elapsed: ' num2str(timeElapsed)]);
%xlabel(num2str(2*width))
%ylabel(num2str(2*width))
axis equal
xlim([floor(middlex - radius),ceil(middlex + radius)])
ylim([floor(middley - radius),ceil(middley + radius)])

%% Display Outputs

n = particleNumber;
disp(['Half number of particles: ' num2str(halfParticleNumber)]);
disp(['Number of particles: ' num2str(particleNumber)]);
disp(['Fractal Dimension: ' num2str(fractdim)]);
disp(['Radius: ' num2str(distanceFromCenter)]);
disp(['Time Elapsed: ' num2str(timeElapsed)]);
disp(['Double Time: ' num2str(doubleTime)]);




