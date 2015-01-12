%% DLA SIMULATION USING RECTANGULAR MATRIX

clear
clf
clc

% Set up variables
xwidth = input('x-width: ');
ywidth = input('y-width: ');

% numberOfParticles = input('Number of Particles: ');

% Initial Setup:
% We create a matrix of zeros and then aggegate the whole bottom row.
% Change entry of matrix to 0.1-1 if it is aggregated.

matrix = zeros(ywidth,xwidth);

for i = 1:xwidth
    matrix(ywidth,i) = 0.1;
end


% We want to let a particle randomly enter from the the top of the matrix

particleNumber = xwidth; %integer
x = 0; %integer
y = 0; %integer
aggregate = 0; %bool
endScript = 0; %bool



while endScript == 0
    particleNumber = particleNumber+1;
    %first decide where the particle starts from
    %start on top
    randPlace = rand*xwidth;
    randPlace = ceil(randPlace);
    x = randPlace;
    y = 1;
    
    % Now we have a starting point on the top edge of the matrix. Now we want
    % to perform random walk in the matrix until we are next to an entry
    % which is not 0. Set aggregate to false.
    
    aggregate = 0;
    
    while aggregate == 0
        randWalk = rand;
        if x == 1
            % top left corner
            if y == 1
                if randWalk < 0.5
                    x = x+1;
                else
                    y = y+1;
                end
                % bottom left corner
            elseif y == ywidth
                if randWalk < 0.5
                    x = x+1;
                else
                    y = y-1;
                end
                % left edge other than corners
            else
                if randWalk < 1/3
                    x = x+1;
                elseif randWalk < 2/3
                    y = y+1;
                else
                    y = y-1;
                end
            end
        elseif x == xwidth
            % top right corner
            if y == 1
                if randWalk < 0.5
                    x = x-1;
                else
                    y = y+1;
                end
                % bottom left corner
            elseif y == ywidth
                if randWalk < 0.5
                    x = x-1;
                else
                    y = y-1;
                end
                % right edge other than corners
            else
                if randWalk < 1/3
                    x = x-1;
                elseif randWalk <2/3
                    y = y+1;
                else
                    y = y-1;
                end
            end
            % not on left or right edges
        else
            % on top edge but not corners
            if y == 1
                if randWalk < 1/3
                    x = x+1;
                elseif randWalk < 2/3
                    x = x-1;
                else
                    y = y+1;
                end
                % on bottom edge but not corners
            elseif y == ywidth
                if randWalk < 1/3
                    x = x+1;
                elseif randWalk < 2/3
                    x = x-1;
                else
                    y = y-1;
                end
                % just in the middle somewhere!
            else
                if randWalk < 0.25
                    x = x+1;
                elseif randWalk < 0.5
                    x = x-1;
                elseif randWalk < 0.75
                    y = y+1;
                else
                    y=y-1;
                end
            end
        end
        
        
        
        
        
        % now check whether it is now adjacent to a 1
        
        if x == 1
            % top left corner
            if y == 1
                if (matrix(y,x+1) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % bottom left corner
            elseif y == ywidth
                if (matrix(y,x+1) + matrix(y-1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % left edge other than corners
            else
                if (matrix(y,x+1) + matrix(y-1,x) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
            end
        elseif x == xwidth
            % top right corner
            if y == 1
                if (matrix(y,x-1) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % bottom left corner
            elseif y == ywidth
                if (matrix(y,x-1) + matrix(y-1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % right edge other than corners
            else
                if (matrix(y,x-1) + matrix(y-1,x) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
            end
            % not on left or right edges
        else
            % on top edge but not corners
            if y == 1
                if (matrix(y,x+1) + matrix(y,x-1) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % on bottom edge but not corners
            elseif y == ywidth
                if (matrix(y,x+1) + matrix(y,x-1) + matrix(y-1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
                % just in the middle somewhere!
            else
                if (matrix(y,x+1) + matrix(y,x-1) + matrix(y-1,x) + matrix(y+1,x))~=0
                    % Then there is an adjacent particle, so we
                    % need to stick it
                    aggregate=1;
                end
            end
        end
        
        
        
        %We keep repeating this until the particle is finally stuck
        
    end
    
    if aggregate
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
        
        % uncomment this if want to watch diagram 'grow'
        %         figure(1)
        %         hold on
        %         imshow(matrix)
        
        
        % If the particle hits the top then we stop simulation
        % note that, at this point, the variables x and y represent a stuck particle already
        if (y == 1)
            endScript = 1;
            % now want to calculate diameter of fractal
            % we could change to make sure that we always add particles until we hit
            % the edge, then diameter we could work out
            %             xdif = abs(x - middlex);
            %             ydif = abs(y - middley);
            %             radius = sqrt(xdif^2 + ydif^2);
            %             diameter = 2*radius;
            %             fractdim = log(particleNumber)/log(radius);
            
            % Set the escape parameter equal to true to stop the loop
        end % end escape check
    end
    
end

%Tell why escaped, if true
% if (endScript == 1)
%     disp ('Script ended early as aggregation was about to exceed matrix dimensions');
%
% end

% uncomment this for diagram to just appear at end

figure(1)
imagesc(matrix)
colormap(jet)
title(['DLA with ' num2str(particleNumber) ' particles'])
xlabel(num2str(xwidth))
ylabel(num2str(ywidth))
axis equal
xlim([0,xwidth])
ylim([0,ywidth])
disp(['Number of particles: ' num2str(particleNumber)])
% disp(['Fractal Dimension: ' num2str(fractdim)])
% disp(['Diameter: ' num2str(diameter)])



