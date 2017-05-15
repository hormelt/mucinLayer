function [ imMom, layerMom ] = getMoment( dataDir, channel, totFrames, momentOrder )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Inputs and Defaults Options

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = pwd;
end

if ~exist('channel','var') || isempty(channel)
    channel = 1;
end

if ~exist('momentOrder','var') || isempty(momentOrder)
    momentOrder = 2;
end

%% Load images

tifstruct = dir([dataDir,'\*.tif']); %Only want .tifs from FIJI output
framepat = 't\d+'; %works for normal FIJI output

%since stacks can have a lot of images, initialize cell arrays for speed

frames = cell(numel(tifstruct),1);

for j = 1:numel(tifstruct)
    frames(j) = regexp(tifstruct(j).name,framepat,'match');
end

if ~exist('totFrames','var') || isempty(totFrames)
    totFrames = numel(unique(frames));
end

dims = loadZSeriesWithChannels( 1, [], channel, 1);

zheight = size(dims,4);

layerMom = zeros(totFrames,zheight);
imMom = zeros(totFrames,1);

%% calcs 

h = waitbar(0,'Determining intensities...');

for j = 1:totFrames

im = loadZSeriesWithChannels( j, [], 1, 1);

Rmat = levelPlane(im); %gets rotation matrix that will place cells on horizontal x-y plane

im = im(:,:,:,:,1);

im = squeeze(im);

[I1, I2, I3] = ind2sub(size(im), 1:numel(im)); %one frame, channel at a time;

inds = [I1;I2;I3];

rot = @(rotmat, I3)(rotmat(3,1)*inds(1,:)+rotmat(3,2)*inds(2,:)+rotmat(3,3)*inds(3,:));

rotz = rot(Rmat,inds); %these are the z heights in the rotated frame

level = graythresh(im); %Threshhold to find pixels that actually contain mucin (Ohtsu)

im(im(:)<level)=0;

b = [im(:)'; rotz];

b(2,:) = round(b(2,:)); %rounding to turn heights into indices over which we will average or sum

b(:,b(2,:)==0)=[];

for k = 1:zheight
    layerMom(j,k) = moment(b(1,b(2,:)==k),momentOrder);
end

imMom(j) = moment(b(1,b(1,:)>0),momentOrder);

waitbar(j / totFrames)

end
close(h)

% dataOut(:,:,1) = layerLum;
% dataOut(:,:,2) = intLum;

end

