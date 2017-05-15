% subtractplane.m
% Routine to subtract a plane from z = f(x, y);  x and y are row and column
% indices.
% Can optionally maintain or get rid of the zeroth order term in the fit 
%   (i.e. subtract the entire plane), or just subtract the slope of the
%    plane, keeping the mean value fixed, if input keepconst==true
% If image is not M x N x 1, average over 3rd (color) dimension
% Optional input:  matrix planecoeff -- use these as the matrix coeff. to 
%    subtract, rather than calling fitplane.m to calculate
% Output: zsub = plane subtracted matrix
% Output: planecoeffout = calculated plane coefficients

% Calls fitplane.m (% finds best fit for z = a1*x + a2*y + c = z)

% Raghuveer Parthasarathy Mar. 2004
% last modified August 16, 2007 (coefficient input & output)

function [zsub planecoeffout] = subtractplane(z, keepconst, planecoeff)


isuint8 = isa(z, 'uint8'); % is z a 'uint8' variable?
isuint16 = isa(z, 'uint16'); % is z a 'uint16' variable?

% make array z class double
if ~isa(z, 'double')
    z = double(z);
end

if (size(z,3) > 1)
    z = mean(z,3); 
    disp('Color image:  averaging over third (color) dimension!');
end

mz = mean(z(:));
% fit a plane to z, or use input coefficients
if (nargin < 3)
    [a1, a2, c] = fitplane(z);  % fit a plane
    planecoeffout = [a1 a2 c];
else
    planecoeffout = planecoeff;
end    
sz = size(z);
jmat = repmat((1:sz(1))',1,sz(2));
kmat = repmat(1:sz(2), sz(1), 1);
zsub = z - planecoeffout(1)*jmat - planecoeffout(2)*kmat - planecoeffout(3);
mzsub = mean(zsub(:));
if keepconst
    zsub = zsub - mzsub + mz;
end

% if the original array z is uint8, make zsub uint8
if (isuint8)
    zsub = uint8(zsub);
elseif (isuint16)
    zsub = uint16(zsub);
end
