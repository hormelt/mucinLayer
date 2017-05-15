function [ layerLum, intLum ] = getABIntensity( dataDir, channel, totFrames, rescale )
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

dims = loadZSeriesWithChannels( 1, [], 1:3, 1);

zheight = size(dims,4);

layerLum = zeros(totFrames,zheight);
intLum = zeros(totFrames,zheight);

%% calcs 

h = waitbar(0,'Determining intensities...');

for j = 1:totFrames

im = loadZSeriesWithChannels( j, [], 1:3, 1)*rescale;

Rmat = levelPlane(im); %gets rotation matrix that will place cells on horizontal x-y plane

im = im(:,:,:,:,1);

im = squeeze(im);

[I1, I2, I3] = ind2sub(size(im), 1:numel(im)); %one frame, channel at a time;

inds = [I1;I2;I3];

rot = @(rotmat, I3)(rotmat(3,1)*inds(1,:)+rotmat(3,2)*inds(2,:)+rotmat(3,3)*inds(3,:));

rotz = rot(Rmat,inds); %these are the z heights in the rotated frame

level = graythresh(im); %Threshhold to find pixels that actually contain mucin (Ohtsu)

imfull = im; %will not threshold for luminosity profiles at each time step (aesthetic choice)

im(im(:)<level)=0;

b = [im(:)'; rotz];

b(2,:) = round(b(2,:)); %rounding to turn heights into indices over which we will average or sum

b(:,b(2,:)==0)=[];

c = [imfull(:)'; rotz];

c(2,:) = round(c(2,:));

c(:,c(2,:)==0)=[];

for k = 1:zheight
    layerLum(j,k) = mean(c(1,c(2,:)==k)); %luminosity profile by timestep
    intLum(j,k) = sum(b(1,b(2,:)==k)); %integrated intensity
end

waitbar(j / totFrames)

end
close(h)

intLum = sum(intLum,2);
% dataOut(:,:,1) = layerLum;
% dataOut(:,:,2) = intLum;

end