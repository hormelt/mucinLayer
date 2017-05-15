function [ rheoarray ] = makeRheologyArray( filename, rheoarray )
%Creates structured arrays for rheology project functions
%
%INPUTS:
%   filename: string array with entries corresponding to path for data
%       files,
%   rheoarray: structured array with same form as output; input this if you
%       want to add to a structured array instead of starting a new one
%
%OUTPUTS:
%   rheo array: structured array containing data for rheology experiments

if ~exist('rheoarray', 'var') || isempty(rheoarray)
    startind = 1;
else
    startind = numel(rheoarray)+1;
end

for j = startind:numel(filename)
    
    load(filename{j})

    if ~exist('objs_out', 'var') || isempty(objs_out)
        rheoarray(j).objs = NaN;
    else
        rheoarray(j).objs = objs_out; %this is dedrifted? or what? change this input to be what you are using to calculate D, if it is not
    end
    
    clear objs_out
    
    if ~exist('vesicle_radius', 'var') || isempty(vesicle_radius)
        rheoarray(j).vesicle_radius = NaN;
    else
        rheoarray(j).vesicle_radius = vesicle_radius;
    end
    
    clear vesicle_radius
    
    if ~exist('tension', 'var') || isempty(tension)
        rheoarray(j).tension = NaN;
    else
        rheoarray(j).tension = tension;
    end
    
    clear tension
    
    if ~exist('date', 'var') || isempty(date)
        rheoarray(j).date = NaN;
    else
        rheoarray(j).date = date; %I think this should be the address of the file
    end
    
    clear date
    
    if ~exist('temperature', 'var') || isempty(temperature)
        rheoarray(j).temperature = NaN;
    else
        rheoarray(j).temperature = temperature;
    end
    
    clear temperature
    
    if ~exist('composition', 'var') || isempty(composition)
        rheoarray(j).composition = NaN;
    else
        rheoarray(j).composition = composition;
    end
    
    clear composition
    
    if ~exist('fluorophore', 'var') || isempty(fluorophore)
        rheoarray(j).fluorophore = NaN;
    else
        rheoarray(j).fluorophore = fluorophore;
    end
    
    clear fluorophore
    
    if ~exist('scale', 'var') || isempty(scale)
        rheoarray(j).scale = NaN;
    else
        rheoarray(j).scale = scale;
    end
    
    clear scale
    
     if ~exist('timestep', 'var') || isempty(timestep)
         if ~exist('fps', 'var') || isempty(fps)
        rheoarray(j).timestep = NaN;
         else
             rheoarray(j).timestep = 1/fps;
         end
    else
        rheoarray(j).timestep = timestep;
     end
    
     clear timestep
     clear fps
    
    if ~exist('centerposition', 'var') || isempty(centerposition)
        rheoarray(j).centerposition = NaN;
    else
        rheoarray(j).centerposition = centerposition;
    end
    
    clear centerposition
    
    if ~exist('objs_corrected', 'var') || isempty(objs_corrected)
        rheoarray(j).objs_corrected = NaN;
    else
        rheoarray(j).objs_corrected = objs_corrected;
    end
    
    
    
    clear objs_corrected

end

