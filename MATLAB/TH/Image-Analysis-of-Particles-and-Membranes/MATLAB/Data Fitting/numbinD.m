% numbinD.m
%
% bin translational and rotational diffusion coefficients, using 
% logarithmically spaced bins
%
% Inputs
%   Dr : array of rotational diffusion coefficients, e.g. from msadtr.m
%        Row 1 = Dr value, Row 2 = uncertainties (std)
%   Dt : array of rotational diffusion coefficients, e.g. from msadtr.m
%        Row 1 = Dt value, Row 2 = uncertainties (std)
%   num : number of bins. Bins are logarithmically spaced.
%   cut : don't include the first cut data points
%   plotopt : if true, make plots (default false)
%
% Outputs
%   meanDr : weighted mean Dr value in each bin (NaN if empty)
%   stdDr : weighted std. dev. of Dr in each bin (NaN if empty)
%   meanDt : weighted mean Dt value in each bin (NaN if empty)
%   stdDt : weighted std. dev. of Dt in each bin (NaN if empty)
%   nx : number of points that fell in each bin
%
% calls wmean.m for weighted averages
% calls errorxy.m for plotting error bars in x and y
% overall structure lifted from RP's binavg.m
%
% Tristan Hormel
% Last modified 10-17-13

function [meanx, stdx, meany, stdy] = numbinD(x, y, num, plotopt)

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end

% get rid of NaN's

goodind = ~(isnan(x(1,:)) | isnan(x(2,:)) | isnan(y(1,:)) | isnan(y(2,:)));

meanxclean = x(1,goodind); stdxclean = x(2,goodind); meanyclean = y(1,goodind); stdyclean = y(2,goodind);

x = [meanxclean; stdxclean]; y = [meanyclean; stdyclean];

% sort by x

[xsort, xind] = sort(x(1,:),'descend'); %binning by x-coordinate
xsort(2,:) = x(2,xind);
ysort(1,:) = y(1,xind);
ysort(2,:) = y(2,xind);

nrbins = floor(length(xsort(1,:))/num);

meanx = zeros(1,nrbins); meany = zeros(1,nrbins); stdx = zeros(1,nrbins); stdy = zeros(1,nrbins);

% weighted average

for j = 1:nrbins-1
    [meanx(j), stdw, std_simple] = wmean(xsort(1,(num*(j-1)+1):(num*j)), xsort(2,(num*(j-1)+1):(num*j)));
    stdx(j) = sqrt(stdw^2 + std_simple^2/num);
%     stdx(j) = std_simple/sqrt(num);
%     meanDr(j) = mean(Drsort(1,(num*(j-1)+1):(num*j)));
    [meany(j), stdw, std_simple] = wmean(ysort(1,(num*(j-1)+1):(num*j)), ysort(2,(num*(j-1)+1):(num*j)));
    stdy(j) = sqrt(stdw^2 + std_simple^2/num);
%     stdy(j) = std_simple/sqrt(num);

%     meanDt(j) = mean(Dtsort(1,(num*(j-1)+1):(num*j)));
end

[meanx(end), ~, std_simple] = wmean(xsort(1,num*(nrbins-1):end), xsort(2,num*(nrbins-1):end));
stdx(end) = std_simple/sqrt(numel(xsort(1,(num*(nrbins-1):end)))); %kinda clunky but this works
[meany(end), ~, std_simple] = wmean(ysort(1,num*(nrbins-1):end), ysort(2,num*(nrbins-1):end));
stdy(end) = std_simple/sqrt(numel(ysort(1,(num*(nrbins-1):end))));

if plotopt
    warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
    figure('name', 'Binned with even numbers');
    xlabel('x')
    ylabel('y')
    errorxy([x(1,:)' y(1,:)' x(2,:)' y(2,:)'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', .7.*[1 1 1], 'ScaleX', 'log', 'ScaleY', 'log');
    hold on
    errorxy([meanx' meany' stdx' stdy'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', [0.7 0.3 0], 'FaceColor', [0.9 0.6 0.3], 'MarkSize', 8, 'WidthEB', 1.5, 'ScaleX', 'log', 'ScaleY', 'log');
    warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');
end
