function [matrix,particleNumber,diameter] = DLALatticeNoiseReduction(numberOfParticles,noiseParameter);
%% DLA SIMULATION USING SQUARE MATRIX FOR ON LATTICE

tic

%% Initial Setup:
% We create a matrix of zeros and then plant a seed at the centre of the
% matrix. Change entry of matrix to 1 if it is aggregated.
%
% This is the matrix set up, we start from inner square which grows as the cluster grows and consider
% particle escaped if it reaches the outer square. We apply noise reduction
% in the following manner: each growth site has a counter which increases by a
% unit every time the site is visited. An empty site is occupied only after
% its selection M times.
%  -----------
% |    ----    |
% |   |    |   |
% |   |    |   |
% |    ----    |
%  -----------

matrix = zeros(round(numberOfParticles/4));
siteCounter = zeros(round(numberOfParticles/4));

middle = round(numberOfParticles/8);

matrix(middle,middle) = 1;

% Initialise our variables

particleNumber = 1; %integer
endScript = 0; %bool
maximumDistance = 0;

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
    % which is 1. Set aggregate & escape to false.
    
    distanceFromCenter = R;
    aggregate = 0; %bool to exit below while loop
    escape = 0; %bool to exit below while loop
    numberOfSteps = 0; %record number of steps taken until aggregated
    
    while (aggregate == 0) && (escape == 0)
        
        randWalk = rand;
        if distanceFromCenter > 4/3*R
            x0 = x - middle;
            y0 = middle - y;
            r0 = sqrt(x0^2 + y0^2);
            V = ((r0 - R)/(r0 + R))*tan(pi*randWalk);
            x = (R/r0)*(((1-V^2)*x0 - 2*V*y0)/(1+V^2));
            y = (R/r0)*(((1-V^2)*y0 + 2*V*x0)/(1+V^2));
            x = round(x + middle);
            y = round(middle - y);
        else
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
        end
        
        xdistanceFromCenter = abs(x - middle);
        ydistanceFromCenter = abs(y - middle);
        distanceFromCenter = sqrt(xdistanceFromCenter^2 + ydistanceFromCenter^2);
        
        
        % aggregate test
        if (matrix(y,x+1) + matrix(y,x-1) + matrix(y+1,x) + matrix(y-1,x)) ~= 0
            aggregate = 1;
        end
        
        
        %We keep repeating this until the particle is finally stuck or
        %escapes
        
    end
    
    %if escape
    %disp('particle escaped')
    %end
    
    if aggregate
        siteCounter(y,x) = siteCounter(y,x) + 1;
        if siteCounter(y,x) == noiseParameter;
            particleNumber = particleNumber + 1;
            stuck_particles(particleNumber,1) = x;
            stuck_particles(particleNumber,2) = y;
            stuck_particles(particleNumber,3) = distanceFromCenter;
            if distanceFromCenter > maximumDistance
                maximumDistance = distanceFromCenter;
            end
            
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
            
            matrix(y,x) = 1;
            
            % end of script check
            
            if particleNumber == numberOfParticles
                endScript = 1;
                % calculate fractal dimension
                
                diameter = 2*distanceFromCenter;
                
            end
        end
    end
    
end

%% Plot graph

%matrix = matrix/particleNumber;
imagesc(matrix)
colormap(gray)
title(['Noise Reduction DLA with ' num2str(noiseParameter) ' as noise parameter and ' num2str(numberOfParticles) ' particles'])
%text(width/3,3*R/2 + R/6,['Fractal Dimension: ' num2str(fractalDimension)]);
%text(R/3,3*R/2,['Radius: ' num2str(maximumDistance)]);
timeElapsed = toc;
%text(R/3,3*R/2 - R/6, ['Time Elapsed: ' num2str(timeElapsed)]);
%xlabel(num2str(2*width))
%ylabel(num2str(2*width))
%axis equal
miny = min(stuck_particles(:,1));
maxy = max(stuck_particles(:,1));
minx = min(stuck_particles(:,2));
maxx = max(stuck_particles(:,2));
xdif = maxx-minx;
ydif = maxy-miny;
maxaxis = max(xdif,ydif)/2;
xlim([middle - maxaxis - 10,middle + maxaxis + 10])
ylim([middle - maxaxis - 10,middle + maxaxis + 10])
axis off

%% Display Outputs

disp(['Number of particles: ' num2str(particleNumber)]);
%disp(['Fractal Dimension: ' num2str(fractalDimension)]);
disp(['Diameter: ' num2str(diameter)]);
disp(['Time Elapsed: ' num2str(timeElapsed)]);
