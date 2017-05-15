% AVI2_TIFF_im.m
%
% Function to convert an AVI movie to TIFF files, or a matrix of frames
% (e.g. for tracking)
% 
% discards colormap data
%
% Option:  Do particle tracking, using im2obj_rp.m and dependent functions.
%   Must first have determined objsize, threshold, e.g. with testthresh
%   Output:  object matrix (objs), as usual -- see particle tracking
%   functions.  If no tracking, objs = [].
%
% WARNING:  Will collapse color layers into one grayscale image for
% image matrix and for tracking
%
% Raghuveer Parthasarathy
% 25 June 2008
% Last modifed June 30, 2009

function [im objs] = AVI2_TIFF_im


im = [];

disp('  ');
disp('WARNING:  Will collapse color layers into one grayscale image ');
disp('          for image matrix and for tracking');
disp('Should ideally check if image is color and preserve this...');


% ----------------------------------------------------------------------

% Load Movie, and get information
disp(' ');
disp('AVI2_TIFF_im.m');
disp('   Suggestion: first change to movie directory, using "cd [directory in single quotes]", or toolbar.');
disp(' ');

% Options
prompt = {'Enter 1 to choose filenames from a dialog box', ...
    'If manual entry, enter filename (assumes present directory): '};
dlg_title = 'Series option'; num_lines= 1;
def     = {'1','movie.avi'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
loadopt = logical(str2double(answer(1)));
mFileName = char(answer(2));

if (loadopt==1)
    % dialog box for filename
    [mFileName,mPathName] = uigetfile('*.*', 'AVI movie to load...'); 
    fs = sprintf('   Path Name: %s', mPathName); disp(fs);
    fs = sprintf('   File Name: %s', mFileName); disp(fs);
else
    mPathName = pwd;
end
Afile = strcat(mPathName, mFileName);
Nf = length(mFileName);  % number of characters in file name
if ~strcmpi(mFileName(Nf-3:Nf),'.avi')
    disp('Error -- file names must end in ".avi"!  Press Control-C."');
    pause
end

% information about the movie file
fileinfo = aviinfo(Afile);
Nframes = fileinfo.NumFrames;

% Load certain frames
f1 = sprintf('First frame (of %d) to load:', Nframes);
f2 = sprintf('Last frame (of %d) to load:', Nframes);
prompt = {f1, f2};
dlg_title = 'Loading option'; num_lines= 1;
def     = {num2str(1),num2str(Nframes)};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
firstframe = str2double(answer(1));
lastframe = str2double(answer(2));

% load AVI
disp('   Loading movie (wait for "done" indication) ...');
mov = aviread(Afile, firstframe:lastframe);
disp('      ...done.');

% ----------------------------------------------------------------------

% Conversion
prompt = {'Save TIFFs?  (0==no, 1==yes)', ...
    '  If yes, force GRAYSCALE TIFFs?  (0==no, 1==yes)', ...
    'Output image matrix?  (0==no, 1==yes)'};
dlg_title = 'Rotation parameters'; num_lines= 1;
def     = {'0', '1', '0'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
TIFFout = logical(str2double(answer(1)));
grayTIFFout = logical(str2double(answer(2)));
imout = logical(str2double(answer(3)));

if TIFFout
    % for output files
    digits = ceil(log10(Nframes+1));
    formatstr = strcat('%0', num2str(digits), 'd');
end
if imout
    % allocate memory
    im = zeros(fileinfo.Height, fileinfo.Width, lastframe-firstframe+1);
end

% -----------------------------------------------------------------------
% Options for particle tracking
prompt = {'Track objects (w/ im2obj_rp etc) (1==yes):',...
    'if yes: object size (px)', ...
    'if yes: threshold'};
dlg_title = 'Particle tracking'; num_lines= 1;
def     = {'0', '8', '0.997'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
trackopt = logical(str2double(answer(1)));
objsize = str2double(answer(2));
threshold = str2double(answer(3));

% Can't preallocate memory for tracking -- don't know # particles
objs = [];


% -----------------------------------------------------------------------

progtitle = 'Progress converting frames...';
progbar = waitbar(0, progtitle);  % will display progress
for j=1:(lastframe-firstframe+1),
    tempframe = mov(j).cdata;
    if TIFFout
        % write TIFF
        framestr = sprintf(formatstr, j);
        OutFileName = strcat(mFileName(1:Nf-4), framestr, '.tif');
        if grayTIFFout
            % force gray
            imwrite(mean(tempframe,3), OutFileName, 'tif');
        else
            imwrite(tempframe, OutFileName, 'tif');
        end
    end
    if imout
        im(:,:,j) = mean(tempframe,3); % collapsing color layers, if any
    end
    if trackopt
        % Track data, using function im2obj_rp, fo4_rp, bpfilter
        tempobjs = im2obj_rp(mean(tempframe,3), objsize, threshold);
        tempobjs(5,:) = repmat(j,1,size(tempobjs,2));
        objs = [objs tempobjs];
    end
    waitbar(j/(lastframe-firstframe+1), progbar, progtitle);
end
close(progbar);


