% wmean.m
% weighted average
%
% uncertainty assuming uncorrelated error manifested in stdy
%
% output:
%     meany = weighted mean
%     stdw = uncertainties from input st.dev.
%     std_simple = standard deviation of all the input y values, neglecting
%           weights
%     N = number of points
% Raghuveer Parthasarathy
% August 29, 2008
% last modified Oct. 3, 2008

function [meany, stdw, std_simple, N] = wmean(y, stdy)

if (length(y)==length(stdy))
    stdy2 = stdy.*stdy;
    w = 1./stdy2;
    meany = sum(y.*w) / sum(w);
    stdw = (1/sum(w))*sqrt(sum(w.*w.*stdy.*stdy)); % = 1/sqrt(sum(w))
    std_simple = std(y);
    N = length(y);
else
    meany = NaN;
    stdw = NaN;
    std_simple = NaN;
    N = NaN;
    disp('Error!  Arrays not the same length!');
end

