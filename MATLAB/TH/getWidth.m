function [ layerWidth, thresh ] = getWidth( thresh, dataDir, channel, totFrames, scale )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Inputs and Default Options

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = pwd;
end

if ~exist('channel','var') || isempty(channel)
    channel = 1;
end

if ~exist('rescale','var') || isempty(scale)
    scale = 1;
end

%% Load images

tifstruct = dir([pwd,'\*.tif']); %Only want .tifs from FIJI output
framepat = 't\d+'; %works for normal FIJI output

%since stacks can have a lot of images, initialize cell arrays for speed

frames = cell(numel(tifstruct),1);

for j = 1:numel(tifstruct)
    frames(j) = regexp(tifstruct(j).name,framepat,'match');
end

if ~exist('totFrames','var') || isempty(totFrames)
    totFrames = numel(unique(frames));
end

dims = loadZSeriesWithChannels( 1, [], 1, 1);

zheight = size(dims,4);

layerWidth = zeros(1,totFrames);
%% calcs 

h = waitbar(0,'Determining widths...');

% Tommy's suggestion: get threshhold from first image instead of
% calculating a new one at each image. Captures the fact that intensity
% trace is not just broadening but also getting taller.

if ~exist('thresh','var') || isempty(thresh)

im = loadZSeriesWithChannels(1,[],1,1);

im = im(:,:,:,:,1); %There is a bug here, shouldn't have to do this line. It works though, low priority fix

rowMax = zeros(1,size(im,1));

for j = 1:size(im,1)
    thisRow = mean(im(j,:,:),2);
    rowMax(j) = max(thisRow);
end

thresh = .5*mean(rowMax); % keep this same threshhold for everything. Have to think about if this makes everything relative to the first frame

end

for j = 1:totFrames
   
    im = loadZSeriesWithChannels( j, [], 1, 1);

    im = im(:,:,:,:,1);
    
im = squeeze(im);

widthHere = zeros(1,size(im,1));
% maxHere = zeros(size(im,1),size(im,2));
% ind = zeros(size(im,1),size(im,2));
% RHS = zeros(size(im,1),size(im,2));
% LHS = zeros(size(im,1),size(im,2));
for jj = 1:size(im,1)
        thisRow = mean(im(jj,:,:),2);
    [~, ind] = max(thisRow);
    ind = min(ind); % multimple maxes can happen if image is saturated
    [~,LHS] = min(abs(thisRow(1:ind)-thresh));
    [~,RHS] = min(abs(thisRow(ind:end)-thresh));
    LHS = max(LHS); RHS = min(RHS);
widthHere(jj) = (RHS+(ind-LHS))*scale;
end


layerWidth(j) = mean(widthHere(:));

waitbar(j / totFrames)

end

close(h)
end

