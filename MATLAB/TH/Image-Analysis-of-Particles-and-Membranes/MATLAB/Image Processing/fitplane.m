% fitplane.m
% routine to fit a plane to z = f(x, y)
% Raghuveer Parthasarathy, 2003
% last modified June 3, 2008

function [a1, a2, c] = fitplane(z, x, y)

% finds best fit for z = a1*x + a2*y + c = z
% x and y are the row and column positions of the elements of z.
%   If x and y are not specified, use rows and columns of 2D z array. 

% array z should be of class double
if (isa(z, 'double') ~= 1)
    z = double(z);
end

if (nargin < 2)
    sz = size(z);
    x = repmat((1:sz(1))',1,sz(2));
    y = repmat(1:sz(2), sz(1), 1);
end

% sums
sx = sum(x(:));
sy = sum(y(:));
xx = x.*x;
sxx = sum(xx(:));
yy = y.*y;
syy = sum(yy(:));
xy = x.*y;
sxy = sum(xy(:));
N = size(z,1)*size(z,2);

xz = x.*z;
sxz = sum(xz(:));
yz = y.*z;
syz = sum(yz(:));
sz = sum(z(:));

mat = [sxx sxy sx; sxy syy sy; sx sy N];
dmat = det(mat);

deta1mat = det([sxz syz sz; sxy syy sy; sx sy N]);
deta2mat = det([sxx sxy sx; sxz syz sz; sx sy N]);
detcmat = det([sxx sxy sx; sxy syy sy; sxz syz sz]);

a1 = deta1mat / dmat;
a2 = deta2mat / dmat;
c = detcmat / dmat;

