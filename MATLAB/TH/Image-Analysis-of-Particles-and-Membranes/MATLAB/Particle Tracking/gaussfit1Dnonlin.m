% gaussfit1Dnonlin.m
%
% function to fit a 1D gaussian using non-linear regression
%
% simple 1D version of gaussfit2Dnonlin
% UNLIKE gaussfit2Dnonlin, non-uniform x-positions can be INPUT
% fit to form: z = A*exp(-(x-x0)^2  / (2*sigma^2)) + offset
% Assume offset >= 0
%
% the use of the nonlinear regression is based on a post by John D'Errico at
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/165035
%
%
% Input: 
%    x : 1D array (optional -- default = 1:length(z))
%    z : 1D array
%    tolerance in z for fitting (default 1e-3 if empty)
%    params0 : [optional] starting values of the parameters
%              If not input or empty, default values
%              1 - constant offset (minimal value of z)
%              2 - center (the "center of mass" )
%              3 - sigma   (second moment)
%              4 - amplitude (max. z - min. z)
%    LB     : [optional] lower bounds of search parameters
%              If not input or empty, default values
%              1 - constant offset (0)
%              2 - min. value of x
%              3 - sigma   (0)
%              4 - amplitude (0)
%    UB     : [optional] upper bounds of search parameters
%              If not input or empty, default values
%              1 - constant offset (max value of z)
%              2 - max. value of x.
%              3 - sigma   (Inf.)
%              4 - amplitude (Inf.)
%    lsqoptions : [optional] options structure for nonlinear least-squares
%              fitting, from previosly running 
%              "lsqoptions = optimset('lsqnonlin');"
%              Inputting this speeds up the function.
% Outputs
%    A  : Gaussian amplitude
%    x0 : Gaussian center 
%    sigma: Std dev. of Gaussian 
%    offset : constant offset
%
% Raghuveer Parthasarathy 
% December 2, 2011
% Last modified December 2, 2011

function [A, x0, sigma, offset] = gaussfit1Dnonlin(x, z, tolz, params0, LB, UB, lsqoptions)

N = length(z);

% defaults for initial parameter values, and lower and upperbounds
if ~exist('x', 'var') || isempty(x)
    x = 1:N;
end
if ~exist('tolz', 'var') || isempty(tolz)
    tolz = 1e-3;
end
if ~exist('params0', 'var') || isempty(params0)
    % Default center: "center of mass"
    params0_2 = sum(x.*z)/sum(z);
    % Default width: second moment
    x2 = (x-params0_2).*(x-params0_2);
    params0_3 = sqrt(sum(x2.*z)/sum(z))/2;
    params0 = [min(z(:)), params0_2, params0_3, max(z(:))-min(z(:))];
end
if ~exist('LB', 'var') || isempty(LB)
    LB = [0,min(x),0,0];
end
if ~exist('UB', 'var') || isempty(UB)
    UB = [max(z(:)),max(x),Inf,Inf];
end
if ~exist('lsqoptions', 'var') || isempty(lsqoptions)
    lsqoptions = optimset('lsqnonlin');
end

% More fitting options
lsqoptions.TolFun = tolz;  %  % MATLAB default is 1e-6
lsqoptions.TolX = 1e-5';  % default is 1e-6
lsqoptions.Display = 'off'; % 'off' or 'final'; 'iter' for display at each iteration
params = lsqnonlin(@(P) objfun(P,x,z),params0,LB,UB,lsqoptions);
A = params(4);
x0 = params(2);
sigma = params(3);
offset = params(1);

end

    function resids = objfun(params,x,z)
        temp = x(:) - params(2);
        pred = params(1) + params(4)*exp(-sum(temp.*temp,2)/2/params(3)/params(3));
        resids = pred - z(:);
    end

