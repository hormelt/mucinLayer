% make_ellipse.m
% 
% function to create a filled-in ellipse of a given size, position, 
% and orientation
%
% Inputs
% N : image size (N x N pixels)
% ra, rb : semimajor, semiminor axes
% theta : orientation angle relative to x axis, radians
% xc, yc : center, px (default 0.0), relative to image center.  Note that y
%     increases downward
% 
% Outputs
% im : image; 1 inside ellipse; 0 outside
%
% Raghuveer Parthasarathy
% August 9, 2012


function im = make_ellipse(N, ra, rb, theta, xc, yc)

if ~exist('xc', 'var') || isempty(xc)
    xc = 0.0;
end
if ~exist('yc', 'var') || isempty(yc)
    yc = 0.0;
end

im = zeros(N);

% array of positions, relative to center
xp = repmat((1:N)-(N+1)/2-xc, N, 1);
yp = repmat((1:N)'-(N+1)/2-yc, 1, N);
rp = sqrt(xp.*xp + yp.*yp);
% polar angle
phi = atan2(yp, xp);

% equation of an ellipse, polar coordinates
r_ellipse = ra*rb ./ sqrt((rb*cos(phi-theta)).^2 + (ra*sin(phi-theta)).^2);

im(rp <= r_ellipse) = 1;

