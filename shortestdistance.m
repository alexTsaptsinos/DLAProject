function [d,p] = shortestdistance(point,line1,line2)
%% Function to calculate shortest distance between line and point
% Arguments: 
% point = vector of point
% line1 = one endpoint of the line
% line2 = the other endpoint of the line

% Outputs:
% d = shortest distance
% p = point from shortest distance
%

% get equation of line in form ax + by + c = 0

m = (line2(2) - line1(2))/(line2(1) - line1(1));
a = m;
b = -1;
c = line1(2) - m*line1(1);

d = abs(a*point(1) + b*point(2) + c)/sqrt(a^2 + b^2);
x_1 = (b*(b*point(1) - a*point(2))-a*c)/(a^2 + b^2);
y_1 = (a*(-b*point(1) + a*point(2)) - b*c)/(a^2+b^2);
p = [x_1,y_1];






