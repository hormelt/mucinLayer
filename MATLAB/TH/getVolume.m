function [ vol ] = getVolume( dataDir, channel, method, processop, zrange )
%Function to determine the volume of mucin from confocal images using
%thresholding
%   INPUTS-
%dataDir: Directory containing .tif microscope images
%channel: flourescence channel of interest
%method: 2 options, yz or xy, that determine which plane to use for image
%   processing. Default: yz.
%processop: option for bandpass filtering ('bp'). Seems to work better on
%   unprocessed images ('none').
%   zrange: limit z values to consider, if desired.
%
%   OUTPUTS-
%vol- volume of mucin layer (boxels).
%
%   DEPENDENCIES-
%quickBPF.m, imthresh.m

%% Inputs and Defaults Options

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = pwd;
end

if ~exist('channel','var') || isempty(channel)
    channel = 1;
end

if ~exist('method','var') || isempty(method)
    method = 'yz';
end

%% Load images

tifstruct = dir([pwd,'\*.tif']); %Only want .tifs from FIJI output
framepat = 't\d+'; %works for normal FIJI output

%since stacks can have a lot of images, initialize cell arrays for speed

frames = cell(numel(tifstruct),1);

for j = 1:numel(tifstruct)
    frames(j) = regexp(tifstruct(j).name,framepat,'match');
end

totFrames = numel(unique(frames));

%% Calculations

vol = zeros(1,totFrames); %initiliaze variable

for j = 1:totFrames
    
    if ~exist('zrange','var') || isempty(zrange)
        im = loadZSeriesWithChannels( j, [], channel );
    else
        im = loadZSeriesWithChannels( j, zrange, channel );
    end
    
    
    switch lower(method)
        case 'yz'
            stack = squeeze(permute(im,[4,2,3,1,5]));
        case 'xy'
            stack = im;
    end
    
    switch lower(processop)
        case 'bp'
            stack = quickBPF(stack,3,20,[]);
            threshim = imthresh(stack,[],[]);
        case 'none'
            threshim = imthresh(stack,[],[]);
    end
    
    vol(j) = sum(threshim(:)); %count number of pixels in mucin layer
    
end

end

