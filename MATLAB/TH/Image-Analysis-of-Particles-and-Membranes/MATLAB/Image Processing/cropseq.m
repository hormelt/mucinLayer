% cropseq.m
%
% Function to crop a sequence of images, e.g. "Aseq" saved by get_series.
%   Use saveseq.m to save the cropped sequence
% Optional input argument: cropping rectangle [xmin ymin width height];
% Also saves dateSeq(), each element of which is a cell array of strings 
%   (if this is unavailable, use [] as the input for dateSeq).
% Raghuveer Parthasarathy
% June 6, 2007


function [outA outdateSeq] = cropseq(Aseq, dateSeq, rect)

Nframes = length(Aseq);

% Options
prompt = {'Analyze every "k-th" image; k = ', ...
    'Crop images?  (0==no, 1==yes):'};
dlg_title = 'Series option'; num_lines= 1;
def     = {'1','1'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
Nskip = str2double(answer(1));
cropopt = logical(str2double(answer(2)));

% Load first image
Fr1 = Aseq(1).cdata;
imclass = class(Fr1);  % Determine the class -- uint8 or uint16
s = size(Fr1);


% Region to Crop
if cropopt
    % Determine cropping region, based on first image
    disp('   Select the region to keep...');
    % For ease of cropping, scale the image
    dFr1 = double(Fr1);
    scaleFr1 = uint8((dFr1 - min(dFr1(:)))*255.0 / (max(dFr1(:)) - min(dFr1(:))));
    if (nargin > 2)
        % cropping rectangle is specified
        cropFr1 = imcrop(scaleFr1, rect);
    else
        % interactive cropping rectangle
        [cropFr1,rect] = imcrop(scaleFr1);
    end
    % preallocate memory
    if strcmp(imclass, 'uint8')
        outA = zeros(size(cropFr1,1), size(cropFr1,2), length(1:Nskip:Nframes), 'uint8');
    elseif strcmp(imclass, 'uint16')
        outA = zeros(size(cropFr1,1), size(cropFr1,2), length(1:Nskip:Nframes), 'uint16');
    else
        disp('Image must be 8 or 16 bit!  Press Control-C');
        pause
    end
else
    rect = [1 1 (s(2)-1) (s(1)-1)];
    % preallocate memory
    if strcmp(imclass, 'uint8')
        outA = zeros(s(1), s(2), length(1:Nskip:Nframes), 'uint8');
    elseif strcmp(imclass, 'uint16')
        outA = zeros(s(1), s(2), length(1:Nskip:Nframes), 'uint16');
    else
        disp('Image must be 8 or 16 bit!  Press Control-C');
        pause
    end
end
fs = sprintf('Cropping rectangle: %d % d %d %d', rect); disp(fs);


% Load frames and crop (if desired)
outk=1;
frmin = 1; 
frmax = length(Aseq);
progtitle = 'Progress...'; 
progbar = waitbar(0, progtitle);  % will display progress
for k=frmin:Nskip:frmax,
    tempoutA  = Aseq(k).cdata;
    if cropopt
        outA(:,:,outk) = imcrop(tempoutA,rect);
    else
        outA(:,:,outk) = tempoutA;
    end
    outk = outk+1;
    waitbar((k-frmin)/(frmax-frmin), progbar, progtitle);
end
close(progbar);
Noutframes = outk-1;

outdateSeq = dateSeq;

% Information, including Date Stamps
disp(' ');
fs = sprintf('%d output frames, from input frames %d to %d, increment %d', ...
    Noutframes, frmin, frmax, Nskip); disp(fs);
fs = sprintf('Initial Datestamp: %s', char(outdateSeq(1))); disp(fs);
fs = sprintf('Final Datestamp:   %s', char(outdateSeq(Noutframes))); disp(fs);

