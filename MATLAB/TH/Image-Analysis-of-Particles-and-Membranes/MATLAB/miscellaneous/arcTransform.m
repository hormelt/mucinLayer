function [ correctedarray ] = arcTransform( posarray, Cent, Rad, invert)
% Takes x,y positions from a plane and projects them on to the surface of a
% sphere with radius Rad and center Cent.
%INPUTS: 
%   posarray: positions, in format [x; y]. Can just use normal objs array.
%   Cent: Center of sphere, in format [x; y]- if the sphere is moving,
%       input a sphere center corresponding to each x-y position pair.
%   Rad: Sphere radius- can also be the length of x-y pairs, or just 1
%       value.
%   invert: option to invert the transformation (i.e. take coordinates from
%       sphere to plane).
%OUTPUTS:
%   %s: transformed positions on the surface of the sphere, coordinates
%       corresponding to distances (arc lengths) from the sphere pole. Note
%       that this is NOT longitude/latitude coordinates!

%Note: use looping for objs arrays

if ~exist('invert', 'var') || isempty(invert)
    invert = 0; 
end

% meanCent = mean(Cent,2);

pole(1,:) = repmat(Cent(1,1),1,size(posarray,2)); % center of the GUV
pole(2,:) = repmat(Cent(1,2),1,size(posarray,2));

% Note that matlab counts rows starting at the top- this gets accounted for
% when I add pole(2,:) instead of subtracting

transx = (posarray(1,:) - pole(1,:)); % translate coordinates so that that projection point is at the pole
transy = (posarray(2,:) - pole(2,:));

if invert
    
    % take coordinates from sphere and project to flat surface
    
    x = Rad.*cos(posarray(1,:)).*sin(posarray(2,:));
    y = Rad.*sin(posarray(1,:));
% 
%     x = cos(posarray(1,:)).*sin(posarray(2,:));
%     y = sin(posarray(1,:));
    
    correctedarray = [x; y];
    
else

    % from plane (microscope image) to sphere (e.g. GUVs)
    
rho = sqrt(transx.^2 + transy.^2); % convenient quantity

c = asin(rho./Rad); % convenient quantity

% c = asin(rho);

phi = asin(transy.*sin(c)./rho); % latitude
lambda = atan2(rho.*sin(transx.*sin(c)./rho./cos(c)),rho.*cos(transx.*sin(c)./rho./cos(c))); %longitude
% 
% lambda = atan(transx.*sin(c)./rho./cos(c));

correctedarray = posarray;

correctedarray(1,:) = phi; correctedarray(2,:) = lambda;

end


