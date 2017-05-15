function [ cells ] = IDCells( dataDir, featSize, frameRange, zRange, threshcut, useLog, plotopt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Inputs and Defaults Options

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = pwd;
end

if ~exist('channel','var') || isempty(channel)
    channel = 0;
end

%% Load images and initialize variables

tifstruct = dir([pwd,'\*.tif']); %Only want .tifs from NIS, FIJI output
framepat = 't\d+'; %works for normal FIJI, NIS output

%since stacks can have a lot of images, initialize cell arrays for speed

frames = cell(numel(tifstruct),1);

for j = 1:numel(tifstruct)
    frames(j) = regexp(tifstruct(j).name,framepat,'match');
end

totFrames = numel(unique(frames));

if ~exist('frameRange','var') || isempty(frameRange)
    frameRange = 1:totFrames;
end

% Using the first frame to normalize for photobleaching and initialize some
% stuff

data = squeeze(loadZSeriesWithChannels( 1, [], channel ));
mean0 = mean(data(:)); std0 = std(data(:));
maxx = size(data,1); maxy = size(data,2); maxz = size(data,3);

if ~exist('zRange','var') || isempty(zRange)
    minz = 1;
else
    maxz = numel(zRange); minz = min(zRange);
end

nhoods = zeros(maxx,maxy,maxz);
cells = zeros(5,1);


%% Calculations

% General strategy: make use of the fact that we know cells are roughly
% spheres of a certain size. Threshold, then find the brightest point within
% a cell-sized region.

% We will use image dilation to get just a single identification/cell.
% Here, I'm using the dumbest possible shape to dilate to.

[x,y,z] = ndgrid(-featSize:featSize);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=featSize+1); %simple sphere for structuring element

for j = frameRange
    
    if ~exist('zRange','var') || isempty(zRange)
        data = squeeze(loadZSeriesWithChannels( j, [], channel ));
    else
        data = squeeze(loadZSeriesWithChannels( j, zRange, channel ));
    end
    
    data = (data - mean(data(:)))/std(data(:))*std0 + mean0;
    
    b = bpass3D_TA(data,1,featSize+1);
    
    if useLog
        b = log(b);
        b = b/max(b(:))*255;
    end
    
    imdil = imdilate(b,se);
    threshed = b;
    threshed(threshed(:)<(threshcut*max(threshed(:))))=eps;
    
    nhoods(threshed(:)==imdil(:))=1;
    
    
    for step = (minz-1+featSize):(maxz-featSize)
        
        [y, x] = find(nhoods(:,:,step));
        
        if ~isempty(x)
            
            cells = [cells(1,:) x'; cells(2,:) y'; cells(3,:) step*ones(1,length(x)); cells(4,:) j*ones(1,length(x)); cells(5,:) 1:length(x)];
            
        end
        
    end
    
    
    %% Diagnostics
    
    if plotopt
        
        if mod(j,5)==0
            figure; imagesc(max(data,[],3)); colormap('gray')
            hold on
            plot(cells(1,cells(4,:)==j),cells(2,cells(4,:)==j),'mo')
            figure; imagesc(max(imdil,[],3)); colormap('gray')
            hold on
            plot(cells(1,cells(4,:)==j),cells(2,cells(4,:)==j),'mo')
        end
        
    end
    
end

cells(:,1) = [];


