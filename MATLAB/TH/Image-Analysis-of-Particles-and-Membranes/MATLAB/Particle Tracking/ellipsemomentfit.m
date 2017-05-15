% ellipsemomentfit.m
% function to fit an ellipsoid to a particle (2-bead or "rod") image, 
% using centroid finding and moment analysis to determine the orientation,
% etc. -- i.e. treating it as a 2D probability distribution.  
% See e.g. 
%   http://en.wikipedia.org/wiki/Correlation
%   http://en.wikipedia.org/wiki/Bivariate_normal
%
% Input: 
%    z : 2D image
% Outputs:
%    xc, yc : centroid (px, relative to top left corner)
%    eparams.theta : orientation, radians; note y increases downward, so "upside
%            down"
%    eparams.eig1, eig2 : major and minor axis lengths
%    eparams.ecc : eccentricity
%
% Raghuveer Parthasarathy
% based on quickellipsefit.m, Sept. 13, 2011
% March 15, 2012

function [xc yc eparams] = ellipsemomentfit(z)

[ny nx] = size(z);
[px py] = meshgrid(1:nx, 1:ny);  % a grid of coordinates; note that y increases downward

z = double(z); % make process-able
sumz = sum(z(:));

% make all column vectors
px = px(:); py = py(:); z = z(:);

% centroid
xc = sum(z.*px)/sumz;
yc = sum(z.*py)/sumz;

% Variance and covariance
xvar  = sum(z.*(px-xc).*(px-xc))/sumz;
yvar  = sum(z.*(py-yc).*(py-yc))/sumz;
xyvar = sum(z.*(px-xc).*(py-yc))/sumz;
% covar = [xvar xyvar; xyvar yvar];  % covariance matrix
% [V,D] = eig(covar);

% Calculate eigenvalues of the variance-covariance matrix
% (These are the major and minor axes of the best-fit ellipse)
D = sqrt((xvar-yvar).*(xvar-yvar) + 4*xyvar*xyvar);
eparams.eig1 = sqrt(0.5*(xvar+yvar+D));
eparams.eig2 = sqrt(0.5*(xvar+yvar-D));

% Angle w.r.t. x-axis.  Note that y increases downward
eparams.theta = 0.5 * atan(2*xyvar / (xvar-yvar));
% eparams.theta = atan((eparams.eig1-xvar)/xyvar);  % same

% Eccentricity
eparams.ecc = sqrt(1-(eparams.eig2/eparams.eig1)^2);
% Could also use (eig1-eig2)/(eig1+eig2) as a measure of circularity; I
% think this is the "third eccentricity"

