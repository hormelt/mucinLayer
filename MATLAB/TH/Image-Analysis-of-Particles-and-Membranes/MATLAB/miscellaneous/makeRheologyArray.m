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

if ~exist('objs_out', 'var') || isempty(objs_out)
    rheoarray(j).objs = NaN;
else
    rheoarray(j).objs = load(filename,'objs_out'); %this is dedrifted? or what? change this input to be what you are using to calculate D, if it is not
end

if ~exist('vesicle_radius', 'var') || isempty(vesicle_radius)
    rheoarray(j).vesicle_radius = NaN;
else
    rheoarray(j).vesicle_radius = load(filename,'VESICLE RADIUS');
end

if ~exist('tension', 'var') || isempty(tension)
    rheoarray(j).tension = NaN;
else
    rheoarray(j).tension = load(filename,'TENSION');
end

if ~exist('date', 'var') || isempty(date)
    rheoarray(j).date = NaN;
else
    rheoarray(j).date = load(filename,'DATE'); %I think this should be the address of the file
end

if ~exist('temperature', 'var') || isemptry(temperature)
    rheodata(j).temperature = NaN;
else
    rheodata(j).temperature = load(filename,'TEMPERATURE');
end

if ~exist('composition', 'var') || isempty(composition)
    rheodata(j).composition = NaN;
else
    rheodata(j).composition = load(filename,'COMPOSITION');
end

if ~exist('fluorophore', 'var') || isempty(fluorophore)
    rheodata(j).fluorophore = NaN;
else
    rheodata(j).fluorophore = load(filename,'FLUOROPHORE');
end

if ~exist('scale', 'var') || isempty(scale)
    rheodata(j).scale = NaN;
else
    rheodata(j).scale = load(filename,'scale');
end

if ~exist('centerposition', 'var') || isempty(centerposition)
    rheodata(j).centerposition = NaN;
else
    rheodata(j).centerposition = load(filename,'CENTERPOSITION');
end

if ~exist('objs_corrected', 'var') || isempty(objs_corrected)
    rheodata(j).objs_corrected = NaN;
else
    rheodata(j).objs_corrected = load(filename,'objs_corrected');
end

end

