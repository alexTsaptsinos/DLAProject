function [radius, particleNumber, fractalDimension,timeElapsed] = DLAOffLattice(numberOfParticles);

% This is off-lattice DLA simulation. We aim to utilise methods similar to
% Kuijpers, Martin and Ommen in their 2013 paper. We also utilise the
% method of bringing back to the launching radius in one step seen in
% Sanders 2000 paper.
tic

% We begin by setting up some fixed variables for the simulation

minStepSize = 1; % the minimum step size a particle can take
maxInspectedDistance = 5; % the maximum inspected distance
submatrixSize = 2*maxInspectedDistance+1;
particleRadius = 1; % the radius of a particle
hittingDistance = 2*particleRadius + minStepSize + 1;

% Here are variables that will change as the simulation progresses, but we
% set them up with initial values

particleNumber=1;
maximumDistance = 0;

% We will use various matrices to store information

% stuckParticles will store the exact positions of the particles in the
% cluster
stuckParticles = zeros(numberOfParticles,2);

% the clusterGrid will store an on-lattice version of our cluster
clusterGrid = zeros(numberOfParticles);

% the distanceGrid contains the rounded down distances between the cell and
% the nearest occupied cell. Occupied cells are labeled as zero.
distanceGrid = maxInspectedDistance*ones(numberOfParticles);

% the vicinityGrid contains the rounded down distances between the cells
% and the cell in the middle of the grid. We define it at the beginning and
% then never recalculate it.
vicinityGrid = zeros(submatrixSize);

% Let us set up the vicinityGrid

vicinityMiddle = maxInspectedDistance + 1;

for i=1:submatrixSize
    for j=1:2*maxInspectedDistance+1
        xdistance = abs(vicinityMiddle - j);
        ydistance = abs(vicinityMiddle - i);
        distanceFromMiddle = sqrt(xdistance^2 + ydistance^2);
        vicinityGrid(j,i) = min(floor(distanceFromMiddle),maxInspectedDistance);
    end
end

% Let's start the simulation by planting a seed at the centre
middle = round((numberOfParticles+1)/2);
stuckParticles(particleNumber,1) = 0;
stuckParticles(particleNumber,2) = 0; % not needed but just to demonstrate
clusterGrid(middle,middle)=particleNumber;
submatrix = distanceGrid(middle-maxInspectedDistance:middle+maxInspectedDistance,middle-maxInspectedDistance:middle+maxInspectedDistance);
for i = 1:submatrixSize
    for j=1:submatrixSize
        submatrix(j,i) = min(vicinityGrid(j,i),submatrix(j,i));
    end
end
distanceGrid(middle-maxInspectedDistance:middle+maxInspectedDistance,middle-maxInspectedDistance:middle+maxInspectedDistance) = submatrix;


% Now we have planted our seed at the centre of the matrix and our distance
% grid is initially set up. We can now launch our particles from our launching
% radius which is set to 5 more than the radius of the aggregate.

endscript = 0; %bool to determine when to end the script

while endscript == 0
    
    % First we must launch the particle from the launching radius
    
    Rlaunch = maximumDistance + 5;
    randAngle = 2*pi*rand;
    
    if randAngle < pi/2
        x = Rlaunch*sin(randAngle);
        y = Rlaunch*cos(randAngle);
    elseif randAngle < pi
        theta = randAngle - pi/2;
        x = Rlaunch*cos(theta);
        y = -Rlaunch*sin(theta);
    elseif randAngle < 3*pi/2
        theta = randAngle - pi;
        x = -Rlaunch*sin(theta);
        y = -Rlaunch*cos(theta);
    else
        theta = randAngle - 3*pi/2;
        x = -Rlaunch*cos(theta);
        y = Rlaunch*sin(theta);
    end
    
    % We have where the particle starts from so now we can take our first
    % step of minStepSize
    
    randAngle = 2*pi*rand;
    if randAngle < pi/2
        x = x + minStepSize*sin(randAngle);
        y = y + minStepSize*cos(randAngle);
    elseif randAngle < pi
        theta = randAngle - pi/2;
        x = x + minStepSize*cos(theta);
        y = y - minStepSize*sin(theta);
    elseif randAngle < 3*pi/2
        theta = randAngle - pi;
        x = x - minStepSize*sin(theta);
        y = x - minStepSize*cos(theta);
    else
        theta = randAngle - 3*pi/2;
        x = x - minStepSize*cos(theta);
        y = y + minStepSize*sin(theta);
    end
    
    % now we can begin our walk loop
    aggregate = 0; % bool to break below walk loop
    
    while aggregate == 0
        
        distanceFromCenter = sqrt(x^2 + y^2);
        randAngle = 2*pi*rand;
        if distanceFromCenter > Rlaunch
            % particle has gone outside the launching radius but with method
            % from Sander we can bring it back in one step. Once we have
            % done that we take a step of size minStepSize as if we were
            % beginning from Rlaunch again.
            
            randWalk = rand;
            V = ((distanceFromCenter - Rlaunch)/(distanceFromCenter + Rlaunch))*tan(pi*randWalk);
            x = (Rlaunch/distanceFromCenter)*(((1-V^2)*x - 2*V*y)/(1+V^2));
            y = (Rlaunch/distanceFromCenter)*(((1-V^2)*y + 2*V*x)/(1+V^2));
            
            stepSize = minStepSize;
            
        else
            % the walker is not outside the launching radius so we need to
            % check how far we can walk and whether we will collide
            
            nearestDistance = distanceGrid(middle - round(y),middle + round(x));
            
            if (nearestDistance < hittingDistance) || (nearestDistance == hittingDistance)
                % then walker is so close to the nearest particle in the cluster
                % that there is a chance of collision at the next step
                
                xround = middle + round(x);
                yround = middle - round(y);
                collisionCandidateMatrix = clusterGrid(yround - hittingDistance:yround + hittingDistance,xround-hittingDistance:xround + hittingDistance);
                
                % this matrix contains the labels of stuck particles near
                % the walker, so let's get the labels
                collisionCandidates = [];
                for i = 1:2*hittingDistance+1
                    for j = 1:2*hittingDistance+1
                        if collisionCandidateMatrix(j,i) ~= 0
                            collisionCandidates = [collisionCandidates collisionCandidateMatrix(j,i)];
                        end
                    end
                end
                
                % now we have the labels we can extract the exact
                % co-ordinates from stuckParticles. now we have the current position of the walker (x,y) as
                % well as the positions of all the candidates. We need to
                % generate our random vector direction and check each
                % candidate to see if there is a collision
                % HELP HERE - DOES ANGLE NEED TO BE LESS THAN pi/2??
                
                stepSize = minStepSize;
                
                for i = 1:length(collisionCandidates)
                    tempLabel = collisionCandidates(i);
                    xcol = stuckParticles(tempLabel,1);
                    ycol = stuckParticles(tempLabel,2);
                    
                    A = 1;
                    B = 2*((cos(randAngle))*(x - xcol) + (sin(randAngle))*(y - ycol));
                    C = (xcol - x)^2 + (ycol - y)^2 - (2*particleRadius)^2;
                    coefficients = [A,B,C];
                    equationRoots = roots(coefficients);
                    
                    for j = 1:2
                        root = equationRoots(j);
                        if (imag(root) == 0) && (root > 0)
                            if (root < minStepSize) || (root == minStepSize)
                                % then we have a collision with the
                                % particle after a root step, so adjust the
                                % stepSize
                                if root < stepSize
                                    stepSize = root;
                                    aggregate = 1;
                                end
                            end
                        end
                    end
                end
                
            elseif (nearestDistance > hittingDistance) && (nearestDistance < maxInspectedDistance)
                % then the walker is between the nearest distance and nearest
                % distance + 1 units far from the closest particle in the cluster.
                % A step of size nearestDistance - 4 can be taken without the need
                % of checking for collision
                
                stepSize = nearestDistance - hittingDistance;
                
            elseif nearestDistance == maxInspectedDistance
                % the walker is at least maxInspectedDistance units far from the
                % closest particle in the cluster. A step of size
                % max(distanceFromCentre - clusterRadius - 4,maxInspectedDistance -
                % 4) can be taken without the need to check for collision.
                
                stepSize = max(distanceFromCenter - maximumDistance - hittingDistance,maxInspectedDistance - hittingDistance);
                
            end
            
        end
        
        % we can finally take the step with the appropriately calculated
        % stepSize
        
        x = x + stepSize*cos(randAngle);
        y = y + stepSize*sin(randAngle);
        
        
        
    end
    
    if aggregate
        particleNumber = particleNumber +1;
        stuckParticles(particleNumber,1) = x;
        stuckParticles(particleNumber,2) = y;
        ymat = middle - round(y);
        xmat = middle + round(x);
        clusterGrid(ymat,xmat) = particleNumber;
        % now we have updated the stuckParticle matrix and clusterGrid we
        % need to update our distanceGrid with the new distances. So we
        % take a submatrix around our new aggregated particle and compare
        % to our vicinityGrid
        
        submatrix = distanceGrid(ymat - maxInspectedDistance:ymat + maxInspectedDistance,xmat - maxInspectedDistance:xmat + maxInspectedDistance);
        for i = 1:submatrixSize
            for j=1:submatrixSize
                submatrix(j,i) = min(vicinityGrid(j,i),submatrix(j,i));
            end
        end
        distanceGrid(ymat - maxInspectedDistance:ymat + maxInspectedDistance,xmat - maxInspectedDistance:xmat + maxInspectedDistance) = submatrix;
        
        distanceFromCenter = sqrt(x^2 + y^2);
        if distanceFromCenter > maximumDistance
            maximumDistance = distanceFromCenter;
        end
        disp(['particle aggregated: ' num2str(particleNumber)])
        
        % end of script check
        if particleNumber == numberOfParticles
            endscript = 1;
        end
        
        
    end
    
    
    
end

% Now we have aggregated all the particles and are storing the exact
% locations of their centres in stuckParticles. All we have to do now is
% plot the particles.
ang=0:0.01:2*pi;
particleCosine = particleRadius*cos(ang);
particleSine = particleRadius*sin(ang);

for i = 1:numberOfParticles
    
    xcent = stuckParticles(i,1);
    ycent = stuckParticles(i,2);
    hold on
    plot(xcent+particleCosine,ycent+particleSine,'r');
    
end

hold off
radius = maximumDistance
particleNumber
fractalDimension = log(particleNumber)/log(distanceFromCenter)
timeElapsed = toc



