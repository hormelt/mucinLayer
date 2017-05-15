function [ layerLum, intLum, percentSignal ] = getLayerIntensity( dataDir, channel, totFrames, rescale )
% Calculates the integrated intensity over each image in a z- and t- series
% Inputs:
%   dataDir: Directory with data
%   channel: fluorescence channel to perform the calculation over
%   totFrames: number of frames over which to perform the calculation
%   rescale: scaling value for image intensity
% Outputs:
%   layerLum: Integrated intensity for each layer at each time, stored as
%             time x height
%   intLum: Total integrated intensity over the entire z-stack at each time

%% Inputs and Defaults Options

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = pwd;
end

if ~exist('channel','var') || isempty(channel)
    channel = 1;
end

if ~exist('rescale','var') || isempty(rescale)
    rescale = 1;
end

% keeping generic hard-coded noise filter for now- small filter so should
% be innocuous
lnoise = 2;
kernel_size = 5;
h = fspecial('gaussian', kernel_size, lnoise);

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

layerLum = zeros(totFrames,zheight);
intLum = zeros(totFrames,zheight);

%% Get rotation coordinates

imforrot = loadZSeriesWithChannels( totFrames, [], 1, 1)*rescale; %rotate to final frame (should have most signal)
filtim = double(imfilter(im,h));

imrotatex = squeeze(mean(filtim,1));

xinds = (1:size(filtim,1))';
yinds = (1:size(filtim,1))';
zinds = (1:size(filtim,3));

[xcoord, ycoord, zcoord] = meshgrid(xinds, yinds, zinds);

CoMx = sum(imrotatex.*repmat(zinds,size(imrotatex,1),1),2)./sum(imrotatex,2);

coeffx = pca([xinds CoMx]);

xpcoord = xcoord;
ypcoord = ycoord*coeffx(1,1)-zcoord*coeffx(1,2);
zpcoord = ycoord*coeffx(2,1)+zcoord*coeffx(2,2);

rxim = interp3(xcoord, ycoord, zcoord, filtim, xpcoord, ypcoord, zpcoord);

imrotatey = squeeze(nanmean(rxim,2));

imrotatey(sum(~isnan(imrotatey),2)==0,:,:) = [];
imrotatey(:,sum(~isnan(imrotatey),1)==0,:) = [];

ypinds = ycoord(1:size(imrotatey,1),1,1);
zpinds = (squeeze(zcoord(1,1,1:size(imrotatey,2))))';

CoMy = sum(imrotatey.*repmat(zpinds,size(imrotatey,1),1),2)./sum(imrotatey,2);

coeffy = pca([ypinds CoMy]);


%% calcs 

h = waitbar(0,'Determining intensities...');

for j = 1:totFrames

im = loadZSeriesWithChannels( j, [], 1, 1)*rescale;

im = im(:,:,:,:,1);

im = squeeze(im);
  
% This rotates the imagestack

        xinds = (1:size(filtim,1))';
        yinds = (1:size(filtim,1))';
        zinds = (1:size(filtim,3));
        
        [xcoord, ycoord, zcoord] = meshgrid(xinds, yinds, zinds);
        
        xpcoord = xcoord;
        ypcoord = ycoord*coeffx(1,1)-zcoord*coeffx(1,2);
        zpcoord = ycoord*coeffx(2,1)+zcoord*coeffx(2,2);
        
        rxim = interp3(xcoord, ycoord, zcoord, im, xpcoord, ypcoord, zpcoord);
        
        xppcoord = xcoord*coeffy(1,1)-zcoord*coeffy(1,2);
        yppcoord = ycoord;
        zppcoord = xcoord*coeffy(2,1)+zcoord*coeffy(2,2);
        
        rxyim = interp3(xcoord, ycoord, zcoord, rxim, xppcoord, yppcoord, zppcoord);
        
        rxyim = rxyim(100:(end-30),100:(end-30),:);

level = graythresh(rxyim); %Threshhold to find pixels that actually contain mucin (Ohtsu)

imfull = im; %will not threshold for luminosity profiles at each time step (aesthetic choice)

rxyim(rxyim(:)<level)=0;

percentSignal(j) = sum(rxyim(:)==0)/numel(rxyim);

for k = 1:zheight
    layerLum(j,k) = mean(mean(rxyim(:,:,k))); %luminosity profile by timestep
    intLum(j,k) = sum(sum(rxyim(:,:,k))); %integrated intensity
end

waitbar(j / totFrames)

end
close(h)

intLum = sum(intLum,2);
% dataOut(:,:,1) = layerLum;
% dataOut(:,:,2) = intLum;

end

