% calc_avgprofile.m
%
% subroutine which calculates intensity profiles along lines parallel to
% the central line, and averages them.
% for use with linescan.m
%
% Raghuveer Parthasarathy
% 11 Mar. 04
% last modified Jan. 18, 2011 (remove line drawing)

function [avgc] = calc_avgprofile(A, cx, cy, n, navg)

% returns avgc, the average intensity along cx,cy (of length n) -- double precision
% navg is the number of pixels on each side over which to average
% imgfig is the number of the figure with the image in it, for the line
% plotting.

thetaperp = atan2(-1.0*(cx(n)-cx(1)), cy(n)-cy(1));  % angle of the perpendicular to the linescan line
avgc = zeros(size(cx));
for j=-navg:navg,
    xshift = cx + j*cos(thetaperp);
    yshift = cy + j*sin(thetaperp);
    cpar = improfile(A, xshift, yshift, n);
    avgc = avgc + cpar;
end
avgc = double(avgc) / (2*navg+1);


