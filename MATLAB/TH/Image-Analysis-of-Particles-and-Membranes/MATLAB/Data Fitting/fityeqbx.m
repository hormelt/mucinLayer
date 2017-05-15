% fityeqbx.m
%
% program to fit the equation y = bx to data sets [x], [y] (i.e. a line with
% intercept at the origin)
% incorporate uncertainty in y, or estimate by deviation from fit if sigy = []
%
% Raghu  18 April, 2004
% last modified 10/8/12 : include R^2 values -th

function [B, sigB, chi2, R] = fityeqbx(x, y, sigy)

if (length(x) ~= length(y))
    disp('Error!  x, y are not the same size!')
    qwer = input('Recommend Control-C. [Enter]');
end

N = max(size(x));
sxx = sum(x.*x);
sxy = sum(x.*y);
B = sxy / sxx;
if (nargin < 3)
    sigy = (y - B*x);  % estimate uncertainty by deviation from the line
end
sigB = sqrt(sum(x.*x.*sigy.*sigy))/sxx;
chi2 = sum(sigy.*sigy)/N;
R = 1 - sum((y-B.*x).^2)/sum((y-mean(y)).^2);

