% calc_image_gradient.m
%
% Calculate the gradient magnitude of images, using a Sobel
% filter and locally smoothing.
% Image can be 2D or 3D.  If 3D, can separately examine each 2D slice, or
% do 3D gradient calculation.
% called by markersegment.m
%
% Inputs:
%   A : 2D image, or 3D stack of images.
%   NSobel_and_scaleratio : Sobel filter size; can be 3 (default), 5, or 7
%            If a single number, use this to calculate gradient magnitude
%            for each slice of A separately.
%            If an array of two numbers, do 3D gradient magnitude calcultion.
%            Element #1 is the filter size (presently must be 3), and #2 is
%            the scale ratio in dimension 3 to that of dimension 1, 2.  For
%            example, if the image scale is 0.2 um/px in x, y and 1
%            um/slice in z, then the scale ratio is 5, and NSobel = [3 5]
%   filtwidth : square filter size (for smoothing each gradient component)
%            default 3
%            1 (or smaller) = no filtering
%            For 3D images, if scaleratio>1, do filtering of each 2D slice
%            separately rather than doing a 3D average, since resolution in
%            z is worse.
%
% Raghuveer Parthasarathy
% April 15, 2012 (extracted from markersegment.m)
% March 13, 2013 (3D calculations)
% last modified May 27, 2013


function [gradmag] = calc_image_gradient(A, NSobel_and_scaleratio, filtwidth)

A = double(A);
gradmag = zeros(size(A));
Nframes = size(A, 3);  % number of slices

if ~exist('NSobel_and_scaleratio', 'var') || isempty(NSobel_and_scaleratio)
    NSobel_and_scaleratio = 3;
end
if length(NSobel_and_scaleratio)>1
    scaleratio = NSobel_and_scaleratio(2);
    NSobel = NSobel_and_scaleratio(1);
    calc3Dgradient = true;
else
    % don't do 3D gradient calculation
    calc3Dgradient = false;
    NSobel = NSobel_and_scaleratio;
end

if ~exist('filtwidth', 'var') || isempty(filtwidth)
    filtwidth = 3;
end

if calc3Dgradient && NSobel~=3
    disp('NSobel must be 3 for 3D gradient calculation in calc_image_gradient.m !')
    pause(1)
    disp('Setting to 3');
    NSobel = 3;
end

if filtwidth>1 && mod(filtwidth,2)==0;
        filtwidth = filtwidth+1;  % make sure filter size is odd, if filtering
end

if calc3Dgradient
    % calculate gradient magnitude of 3D image.
    % Presently very clunky, inefficient -- should separate variables,
    % allow different NSobel, etc.!
    
    % 3D 3x3x3 Sobel filter (very clunky, inelegant!), for gradient magnitude.
    hx = zeros([3 3 3]);
    hy = zeros([3 3 3]);
    hz = zeros([3 3 3]);
    hz(:,:,1) = [1 2 1; 2 4 2; 1 2 1];
    hz(:,:,2) = [0 0 0; 0 0 0; 0 0 0];
    hz(:,:,3) = -[1 2 1; 2 4 2; 1 2 1];
    hy(:,:,1) = [1 2 1; 0 0 0; -1 -2 -1];
    hy(:,:,2) = [2 4 2; 0 0 0; -2 -4 -2];
    hy(:,:,3) = [1 2 1; 0 0 0; -1 -2 -1];
    hx(:,:,1) = [-1 0 1; -2 0 2; -1 0 1];
    hx(:,:,2) = [-2 0 2; -4 0 4; -2 0 2];
    hx(:,:,3) = [-1 0 1; -2 0 2; -1 0 1];
    Ax = imfilter(A, hx);
    Ay = imfilter(A, hy);
    Az = imfilter(A, hz)/ scaleratio; % correct for aspect ratio
    % locally average (smooth) the gradient
    if filtwidth>1
        if scaleratio > 1
            % filter each slice with a 2D square filter
            hf = ones(filtwidth, filtwidth)/(filtwidth^2);
            avgAx = imfilter(Ax,hf);  % filter with a square filter
            avgAy = imfilter(Ay,hf);
            gradmag = sqrt(avgAx.^2 + avgAy.^2 + Az.^2);
        else
            % a 3D averaging
            hf = ones(filtwidth, filtwidth, filtwidth)/(filtwidth^3);
            avgAx = imfilter(Ax,hf);  % filter with a square filter
            avgAy = imfilter(Ay,hf);
            avgAz = imfilter(Az,hf);
            gradmag = sqrt(avgAx.^2 + avgAy.^2 + avgAz.^2);
        end
        else
            % no filtering
        gradmag = sqrt(Ax.^2 + Ay.^2 + Az.^2);
    end
else
    % 2D gradient calculation, slice by slice
    hy = sobelN(NSobel);
    hx = hy';
    for j=1:Nframes
        
        % When filtering, note that the input image has to be double precision, not int!
        Ay = imfilter(A(:,:,j), hy, 'replicate');
        Ax = imfilter(A(:,:,j), hx, 'replicate');
        
        % locally average (smooth) the gradient
        if filtwidth>0
            hf = ones(filtwidth)/filtwidth/filtwidth;
            avgAx = filter2(hf,Ax,'same');  % filter with a square filter
            avgAy = filter2(hf,Ay,'same');
            gradmag(:,:,j) = sqrt(avgAx.^2 + avgAy.^2);
        else
            % no filtering
            gradmag(:,:,j) = sqrt(Ax.^2 + Ay.^2);
        end
    end
end

