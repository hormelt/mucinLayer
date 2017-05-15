% culldimobjs.m
%
% function to remove dim objects from an (unlinked) object matrix, 
% e.g. the output of im2obj_rp.m, which calls fo4_rp.m.
% Dim objects can be 
%   (1) brightness less than some integrated brightness threshold, ot
%   (2) dimmer than the n brightest objects in a frame
% Also renumbers objects so that within each frame (row 5 of the object
%   matrix), particle IDs (row 4) are 1, 2, 3, ...
%
% Inputs
%   objs : object matrix.  Note that row 3 is the integrated brightness,
%          returned by fo4_rp.m
%   brthresh : brightness threshold (if option 1), or number of bright 
%          objects to keep (if option 2).
%          Brightness threshold must be previously determined,
%          e.g. by a histogram:  figure; hist(objs(3,:),100)
%   threshopt : which option for identifying dim objects (see above)
%          default (if empty): 1
%
% Output
%   objs_out : new object matrix
%
% Raghuveer Parthasarathy
% March 4, 2011
% Last modified: March 22, 2011

function objs_out = culldimobjs(objs, brthresh, threshopt)

if isempty(threshopt)
    threshopt = 1;
end

if threshopt==1
    objs_out = objs(:,objs(3,:) >= brthresh);  % just the columns with bright objects
elseif threshopt==2
    % keep the n brightest objects in each frame
    utrk = unique(objs(5,:));  % all frame numbers
    objs_out = zeros(size(objs));
    startcol = 1;
    for k=1:length(utrk)
        frcols = find(objs(5,:)==utrk(k));  % indices of columns of this frame
        if length(frcols) <= brthresh
            % there aren't more than "n" objects here, so keep all of them
            lastcol = startcol+length(frcols)-1;
            objs_out(:,startcol:lastcol)=objs(:,frcols);
            startcol = startcol+length(frcols <= brthresh);
        else
            objsk = objs(:,frcols);  % columns corresponding to this frame
            [sbr, sbrix] = sort(objsk(3,:), 2, 'descend');
            lastcol = startcol+1;
            objs_out(:,startcol:lastcol)=objsk(:,sbrix(1:2));
            startcol = startcol+2;
        end
    end
    objs_out = objs_out(:,1:lastcol);
else
    disp('Error!  bad threshold option in culldimobjs.m');
end


% Now re-number
fr = unique(objs_out(5,:));  % all frame numbers
for j=1:length(fr)
    objs_out(4,objs_out(5,:)==fr(j)) = 1:sum(objs_out(5,:)==fr(j));
end

function [ output_args ] = Untitled29( input_args )
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here


end

