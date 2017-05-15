% viewseries.m
%
% Function to view a series of images, either input as an array or read 
%    from a series of files or a multi-page TIFF.  
% Provides an interactive GUI
% Calls getnumfilelist.m to get file names, etc.
% For multiple TIFF files, file names must all be the same *length,* 
%     with the sequence number indicated as trailing digits, e.g.
%     seq001.tif, seq002.tif, ... seq012.tif, ...seq123.tif
% Can read 8- or 16-bit TIFF images.  
% Can read grayscale or COLOR (3-layer RGB) images.
%
% Input: 
%    im : optional, 3D array (series of 2D images)
%    pts : optional, structure of x and y positions to superimpose on the
%          images.  Fields pts(k).x and .y are x and y positions in frame k.
%          Field .plotthispt is an array of 1's and 0's that indicate if
%          each point is to be plotted.  (If this field is empty, ignore)
%          Length of pts must be the same as the number of frames.
% Output: none
%
% Raghuveer Parthasarathy
% June 28, 2009 
% July 10, 2011: Call getnumfilelist and allow multipage TIFF input
% August 22, 2012: allow input of an image array
% November 23, 2012: allow input of an array of points to plot
% last modified March 17, 2013


function viewseries(im, pts)
% nested function

disp(' ');
disp('viewseries.m');
if ~exist('im', 'var')
    im = [];
    disp('   Suggestion: first change to images sequence directory.');
    disp(' ');
end

firstdir = pwd;  % present directory

%% -----------------------------------------------------------------------
% Load images: Get File Names etc. from get numfilelist 
%    Can be multipage TIFF

if ~exist('im', 'var') || isempty(im)
    % Get file names
    [fbase, frmin, frmax, formatstr, FileName1, FileName2, PathName1 ext ismultipage] = ...
        getnumfilelist;
    if isempty(frmin)
        % only one file, so min and max are empty
        frmin = 1;
        frmax = 1;
    end
else
    % array of images is already input
    frmin = 1;
    frmax = size(im,3);
end
Nframes = frmax - frmin + 1;


%% -----------------------------------------------------------------------
% if points are provided, make sure the number matches Nframes
if exist('pts', 'var')
    if length(pts) ~= Nframes
        errordlg('pts array is not the same size as Nframes; ignoring.')
        pts = [];
    end
else
    pts = [];
end

%% -----------------------------------------------------------------------
% Load first image
if isempty(im)
    if ismultipage
        A  = imread(strcat(fbase, ext), 1);
    else
        A = imread(FileName1);
    end
else
    A = im(:,:,1);
end
s = size(A);
iscolor = (length(s)==3);  % TRUE if the image is color (3 layers)

% -------------------------------------------------------------
% GUI -- user interface for viewing images

%  Initialize and hide the GUI as it is being constructed.
fGUI = figure('Name','viewseries', 'Menubar','none', ...
    'Visible','off','Position',[100,100,1000,800], 'Color', [0.2 0.4 0.8]);

% Construct the components.
% Create axes, for the images
hsig = uicontrol('Style','text',...
    'String','viewseries.m -- 2009, R. Parthasarathy',...
    'FontWeight', 'bold', 'Units','normalized', ...
    'Position',[0.8,0.060,0.2,0.040]);
% Initialization
hframetextinit = uicontrol('Style','text',...
    'String','Frame no.', 'FontWeight', 'bold', 'Units','normalized', ...
    'Position',[0.585,755/800,0.085,15/800]);
hframetext  = uicontrol('Style','edit',...
    'String','Frame no.', 'Units','normalized',...
    'Position',[0.675,755/800,0.040,15/800], ...
    'Callback',{@chooseframetext_Callback});
hframe = uicontrol('Style','slider', ...
    'Max', frmax, 'Min', frmin, 'Value', frmin, ...
    'SliderStep', [0.5/Nframes 1/Nframes], 'Units','normalized',...
    'Position', [0.720,755/800,0.140,15/800], ...
    'Callback',{@chooseframe_Callback});
    % Annoying that Mac doesn't have a slider arrow; I'll therefore make
    % the "clicking" step 1 image.
hexit  = uicontrol('Style','pushbutton',...
    'String','Exit', 'Units','normalized',...
    'Position',[0.875,765/800,0.115,30/800],'Callback',{@exit_Callback});
% Messages
hmsg    = uicontrol('Style','text',...
    'String','Message: ', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Units','normalized',...
    'Position',[0.050,780/800,0.400,15/800]);


% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([fGUI, hsig, hframetextinit, hframetext, hframe, hexit, hmsg], ...
    'Units','normalized');

% Create variables here
currframe = frmin;  % current frame being viewed

% Move the GUI to the center of the screen.
movegui(fGUI,'center')
% Make the GUI visible.
set(fGUI,'Visible','on');
figure(fGUI);
set(hframetext, 'String', sprintf('%d', currframe));
set(hframe, 'Value', currframe);
% Show the first image
imagesc(A);
if ~iscolor
    colormap('gray');
end
title('Image');
if ~isempty(pts)
    hold on
    plot(pts(1).x, pts(1).y, 'c+');
end


% Callback functions

    function chooseframetext_Callback(source,eventdata)
        % set frame to view (text entry)
        currframe = round(str2double(get(source,'string')));
        set(hframe, 'Value', currframe);
        loadandshow(currframe);
    end

    function chooseframe_Callback(hObject, source,eventdata)
        currframe = round(get(hObject,'Value'));
        % Also update text entry box:
        set(hframetext, 'String', sprintf('%d', currframe));
        loadandshow(currframe);
    end

    function loadandshow(currframe)
        % load and show image
        if isempty(im)
            if ismultipage
                A  = imread(strcat(fbase, ext), currframe);
            else
                framestr = sprintf(formatstr, currframe);
                FileName = strcat(fbase, framestr, ext);
                A  = imread(FileName);  % image
            end
        else
            A = im(:,:,currframe);
        end
        imagesc(A);
        if ~iscolor
            colormap('gray');
        end
        if ~isempty(pts)
            hold on
            if isfield(pts, 'plotthispt')
                ptind = logical(pts(currframe).plotthispt);
                plot(pts(currframe).x, pts(currframe).y, 'rx');
                plot(pts(currframe).x(ptind), pts(currframe).y(ptind), 'yx');
            else
                plot(pts(currframe).x, pts(currframe).y, 'yx');
            end
        end
    end

    function exit_Callback(source,eventdata)
        % Exit
        cd(firstdir);  %Return to original directory
        close(fGUI)
    end

end