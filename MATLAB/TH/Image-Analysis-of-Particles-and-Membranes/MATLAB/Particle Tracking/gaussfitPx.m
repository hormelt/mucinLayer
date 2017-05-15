% gaussfitPx.m
%
% Gaussian fit of the probability distribution of a 1D variable (x)
% fit to form: P(x) ~ A*exp(-(x-x0)^2 / (2*sigma^2))
% Input:  x -- data points
%         xctrs -- Single value (dx) for bins for histogram, or
%                  array of (xmin+dx/2):dx:(xmax-dx/2) for bin centers for histogram
%                  If empty ([]), use default 21 centers
%         threshold -- Consider only data above some probability threshold,
%            defined as a fraction of the max. value of the histogram of x
%            Default (if only 2 arguments or if threshold==[]): threshold = 0.2  
%            Threshold ensures that fit is not dominated by tails of the Gaussian
%         plotopt -- if true, plot histogram and gaussian fit in new window
%
% Test, e.g., with
%      x = [1.5 + 4*randn(10000,1); 7*randn(3000,1)];  
%      mostly std. dev. of 4, some w. std of 7
%      [x0, sigma] = gaussfitPx(x, 0.5, 0.1, 1)
%
% Raghuveer Parthasarathy Jan. 25, 2011
%    See notes: replaces gaussfit.m for Gaussian prob. distributions
% last modified Feb. 4, 2011

function [x0, sigma] = gaussfitPx(x, xctrs, threshold, plotopt)

x = x(:);   % ensure array shape
if isempty(xctrs)
    Nctrs = 21;
    dx = (max(x)-min(x))/Nctrs;
    xctrs = (min(x)+dx/2):dx:(max(x)-dx/2);
elseif length(xctrs)==1
    % user has input a single number, dx;
    dx = xctrs;
    xctrs = (min(x)+dx/2):dx:(max(x)-dx/2);
else
    % user has input the array
    dx = mean(diff(xctrs));
end

if or((nargin < 3),isempty(threshold))
    threshold = 0.2;
end
if (nargin < 4)
    plotopt=false;
end

if or(length(xctrs)==1, isempty(xctrs))
    disp('Warning: only one bin!  Returning simple mean and std. dev.');
    x0 = mean(x);
    sigma = std(x);
else
    y = hist(x, xctrs); % histogram
    noisex = xctrs(y < threshold*max(y));  % x of bins with low probability
    goody = y(y >= threshold*max(y));  % just for plotting purposes
    
    if isempty(noisex)
        % No noisy bins
        goodx = x;
    else
        isnoisy = zeros(length(noisex),length(x));
        for j=1:length(noisex)
            isnoisy(j,:) = abs(x - noisex(j)) < dx/2.0;  % true if close to this "noisy" bin
        end
        sumisnoisy = sum(isnoisy,1);
        goodx = x(sumisnoisy==0);  % only if not close to any noisy bin
    end
    
    sigma = std(goodx);
    x0 = mean(goodx);
    
    if plotopt
        figure
        plot(xctrs,y/sum(goody)/dx,'ko');
        hold on
        plot(xctrs, (1/sqrt(2*pi*sigma*sigma))...
            *exp(-(xctrs-x0).*(xctrs-x0) / (2*sigma*sigma)), 'r-');
        plot(xctrs, ones(size(xctrs))*threshold/sqrt(2*pi*sigma*sigma), 'b:');
        title('GaussfitPx.m'); xlabel('xctrs'); ylabel('y');
    end
end


