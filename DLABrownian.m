%% DLA Off Latice
% In this simulation, we take the same step size every time then perform a
% test after each step to check if the particle has aggregated/escaped. We
% check if the particle is aggregated by measuring the distance between the
% centre of the particle and the centre of every other particle. If any
% distance is less/equal than 2*radius we know that there must be some
% overlap and so the particle is aggregated. We record the stuck particles
% in a matrix (stuck_particles) which records the x position, y position,
% distance from origin and angle fro origin
clear
clc
clf

%% Set up variables
tic
numberOfParticles = 100;
step_size = 1;
stuck_particles = zeros(1,4);
particleRadius = 1;
ang=0:0.01:2*pi;
particleCosine = particleRadius*cos(ang);
particleSine = particleRadius*sin(ang);



%% Plant initial seed
figure(1);
plot(particleCosine,particleSine,'r');
particleNumber = 1; %integer


x = 0; %integer
y = 0; %integer
endScript = 0; %bool

%% Begin Loop

while particleNumber < numberOfParticles
    % first find initial position on entrance circle
    escape = 0; % bool
    aggregate = 0; % bool
    
    % set up variable starting distance depending on how big the cluster is
    maximumDistance = max(stuck_particles(:,3));
    R = maximumDistance*(3/2); %we try 50% further away from the size of the cluster
    randAngle = 2*pi*rand;
    hold on
    
    if randAngle < pi/2
        x = R*sin(randAngle);
        y = R*cos(randAngle);
    elseif randAngle < pi
        theta = randAngle - pi/2;
        x = R*cos(theta);
        y = -R*sin(theta);
    elseif randAngle < 3*pi/2
        theta = randAngle - pi;
        x = -R*sin(theta);
        y = -R*cos(theta);
    else
        theta = randAngle - 3*pi/2;
        x = -R*cos(theta);
        y = R*sin(theta);
    end
    
    % plot circles for reference, inner circle is starting, outer circle is
    % escape
    
%     xp=R*cos(ang);
%     yp=R*sin(ang);
%     double_xp = 2*R*cos(ang);
%     double_yp = 2*R*sin(ang);
%     plot(xp,yp);
%     plot(double_xp,double_yp);
    
    
    %% Begin Brownian Motion
    
    while (aggregate == 0) && (escape == 0)
        i = i+1;
        randMove = 2*pi*rand;
        
        if randAngle < pi/2
            u = x + step_size*sin(randMove);
            v = y + step_size*cos(randMove);
        elseif randAngle < pi
            theta = randMove - pi/2;
            u = x + step_size*cos(theta);
            v = y - step_size*sin(theta);
        elseif randAngle < 3*pi/2
            theta = randMove - pi;
            u = x - step_size*sin(theta);
            v = y - step_size*cos(theta);
        else
            theta = randMove - 3*pi/2;
            u = x - step_size*cos(theta);
            v = y + step_size*sin(theta);
        end
        % (x,y) is start, (u,v) is end
        
        % escape test
        distanceFromCenter = sqrt(u^2 + v^2);
        if distanceFromCenter > 2*R
            escape = 1;
        end
        
        % aggregate test
        if escape == 0
            
            for i = size(stuck_particles,1):-1:1
                xdif = u - stuck_particles(i,1);
                ydif = v - stuck_particles(i,2);
                particleDistance = sqrt(xdif^2 + ydif^2);
                if particleDistance < 2*particleRadius
                    aggregate = 1;
                    
                    % Need to adjust position of particle to avoid superposition
                    % USING BRAGA/RIBEIRO
                    
                    particleHitX = stuck_particles(i,1);
                    particleHitY = stuck_particles(i,2);
                    A = [u-x,v-y];
                    B = [xdif,ydif];
                    theta = 2*dot(B,A);
                    psi = dot(B,B) - 4*particleRadius^2;
                    chi = dot(A,A);
                    alpha_one = (-theta + sqrt(theta^2 - 4*psi*chi))/(2*chi);
                    alpha_two = (-theta - sqrt(theta^2 - 4*psi*chi))/(2*chi);
                    if abs(alpha_one) < abs(alpha_two)
                        alpha = alpha_one;
                    else
                        alpha = alpha_two;
                    end
                    u = particleHitX + B(1) + alpha*A(1);
                    v = particleHitY + B(2) + alpha*A(2);
                    
                    break
                end
            end
            
        end
        
        
        
        x = u;
        y = v;
    end
    
    if aggregate
        disp('particle aggegated')
        particleNumber = particleNumber + 1;
        stuck_particles(particleNumber,1) = x;
        stuck_particles(particleNumber,2) = y;
        stuck_particles(particleNumber,3) = distanceFromCenter;
        if x > 0
            if y > 0
                particleAngle = atan(x/y);
            elseif y < 0
                particleAngle = atan(y/x) + pi/2;
            else
                particleAngle = pi/2;
            end
        elseif x < 0
            if y < 0
                particleAngle = atan(x/y) + pi;
            elseif y > 0
                particleAngle = atan(y/x) + 3*pi/2;
            else
                particleAngle = 3*pi/2;
            end
        else
            if y < 0
                particleAngle = pi;
            else
                particleAngle = 0;
            end
        end
        stuck_particles(particleNumber,4) = particleAngle;
        plot(x+particleCosine,y+particleSine,'r');
        
    end
    
    if escape
        disp('particle escaped')
    end
    
end

maximumDistance = max(stuck_particles(:,3));
fractalDimension = log(numberOfParticles)/log(maximumDistance);
xlim([-2*R,2*R]);
ylim([-2*R,2*R]);
title(['DLA with ' num2str(particleNumber) ' particles'])
axis square
text(R/3,3*R/2 + R/6,['Fractal Dimension: ' num2str(fractalDimension)]);
text(R/3,3*R/2,['Radius: ' num2str(maximumDistance)]);
timeElapsed = toc;
text(R/3,3*R/2 - R/6, ['Time Elapsed: ' num2str(timeElapsed)]);
hold off

%% Display Outputs

disp(['Number of particles: ' num2str(numberOfParticles)])
disp(['Fractal Dimension: ' num2str(fractalDimension)])
disp(['Radius: ' num2str(maximumDistance)])


