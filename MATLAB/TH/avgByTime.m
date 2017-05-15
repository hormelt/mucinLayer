function [ multiSDs ] = avgByTime( GoodTracks, interval)
%Find the mean of an array between specific intervals
%   Detailed explanation goes here

intRange = min(GoodTracks(3,:)):interval:max(GoodTracks(3,:));

for j = 1:numel(intRange)-1
    multiSDs(1,j) = nanmean(GoodTracks(1,intRange(j)<GoodTracks(3,:)<=intRange(j+1)));
    multiSDs(2,j) = nanstd(GoodTracks(1,intRange(j)<GoodTracks(3,:)<=intRange(j+1)));
    multiSDs(3,j) = intRange(j+1)/2;
end


end

