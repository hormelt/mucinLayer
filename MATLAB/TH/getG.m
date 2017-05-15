function [ Garr ] = getG( GoodTracks,T,a )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

kB = 1.38e-23;

for j = 1:length(GoodTracks)
    
    data = GoodTracks(j).MSDs;
    
    msds = data(:,1)*1e-9*1e-9;
    omega = 2*pi./data(:,2);
    
    alpha = diff(log(msds));
    
    G = 2*kB*T/3/pi/a./msds(2:end)./gamma(1+alpha);
    
    Gp = G.*cos(pi*alpha/2);
    Garr(j).Gp = Gp;
    Gpp = G.*sin(pi*alpha/2);
    Garr(j).Gpp = Gpp;
    
    Garr(j).G = G;
    Garr(j).alpha = alpha;
    Garr(j).omega = omega;
    
end

end

