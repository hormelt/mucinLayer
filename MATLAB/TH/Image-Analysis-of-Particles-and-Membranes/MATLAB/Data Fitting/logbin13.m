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

function [meanDr stdDr meanDt stdDt nx] = numbinD(Dr, Dt, plotopt)

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end

[Drsort, Drind] = sort(Dr(1,:));
Drsort(2,:) = Dr(2,Drind);

[Dtsort, Dtind] = sort(Dt(1,:));
Dtsort(2,:) = Dt(2,Drind);

meanDt = Dtsort(1,1); stdDt = Dtsort(2,1);
meanDr = Drsort(1,1); stdDr = Drsort(2,1);

if plotopt
    warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
    figure('name', 'logbinD');
    xlabel('Dr, radians^2/s')
    ylabel('Dt, \mu m^2/s')
    errorxy([Drsort(1,:)' Dtsort(1,:)' Drsort(2,:)' Dtsort(2,:)'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', 0.7*[1 1 1], 'ScaleX', 'log', 'ScaleY', 'log');
    hold on
    errorxy([meanDr' meanDt' stdDr' stdDt'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
        'EdgeColor', [0.7 0.3 0], 'FaceColor', [0.9 0.6 0.3], 'MarkSize', 8, 'WidthEB', 1.5, 'ScaleX', 'log', 'ScaleY', 'log');
    warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');
end
