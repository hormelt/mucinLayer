% logbinD.m
%
% bin translational and rotational diffusion coefficients, using 
% logarithmically spaced bins
%
% Inputs
%   Dr : array of rotational diffusion coefficients, e.g. from msadtr.m
%        Row 1 = Dr value, Row 2 = uncertainties (std)
%   Dt : array of rotational diffusion coefficients, e.g. from msadtr.m
%        Row 1 = Dt value, Row 2 = uncertainties (std)
%   nbins : number of bins.  Creates logarithmically spaced bins from the
%        min to max value of Dr.
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
% Raghuveer Parthasarathy
% Sept. 1, 2012
% last modified: Nov. 27, 2012 (TH)

function [meanDr, stdDr, meanDt, stdDt, nx, binlb] = logbinD(Dr, Dt, nbins, plotopt)

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end

% bins for Dr
% nbins = 10;
binlb = logspace(log10(0.75*min(Dr(1,:))), log10(1.33*max(Dr(1,:))), nbins+1);

ind=zeros(1,size(Dr,2));
for j=1:nbins+1
    ind=ind+(Dr(1,:)>=binlb(j));
end

meanDr   = zeros(1,nbins);
stdDr = zeros(1,nbins);
meanDt   = zeros(1,nbins);
stdDt = zeros(1,nbins);
nx = zeros(1,nbins);
for k = 1:nbins,
    Drk = Dr(1,ind==k);
    Dtk = Dt(1,ind==k);
    nx(k) = sum(ind==k);
if and(~isempty(Drk),sum(~isnan(Drk)>0))  % not empty, and not all NaN
        % weighted mean and standard deviation. Note that taking a simple
        % average of the std underestimates the error
        [meanDr(k) stdw std_simple] = wmean(Drk, Dr(2,ind==k));
%         meanDr(k) = mean(Drk);
        stdDr(k) = std_simple/sqrt(length(Drk)); %weighted std added in quadrature with sdom
        [meanDt(k) stdw std_simple] = wmean(Dtk, Dt(2,ind==k));
%         meanDt(k) = mean(Dtk);
        stdDt(k) = std_simple/sqrt(length(Drk));
    else
        meanDr(k) = NaN;
        meanDt(k) = NaN;
        stdDr(k) = NaN;
        stdDt(k) = NaN;
    end
end


if plotopt
    warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
    figure('name', 'logbinD');
    xlabel('Dr, radians^2/s')
    ylabel('Dt, \mu m^2/s')
    errorxy([Dr(1,:)' Dt(1,:)' Dr(2,:)' Dt(2,:)'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', 0.7*[1 1 1], 'ScaleX', 'log', 'ScaleY', 'log');
    hold on
    errorxy([meanDr' meanDt' stdDr' stdDt'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', [0.7 0.3 0], 'FaceColor', [0.9 0.6 0.3], 'MarkSize', 8, 'WidthEB', 1.5, 'ScaleX', 'log', 'ScaleY', 'log');
    warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');
end
