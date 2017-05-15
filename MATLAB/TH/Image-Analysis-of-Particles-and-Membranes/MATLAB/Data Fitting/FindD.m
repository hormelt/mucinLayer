function [ Dt, Dr, chi2 ] = FindD( objs, timestep, scale, method, geometry, rotating, verificationIDs, plotopt, Rad, mindt, cutoff)
%Estimates Diffusion Coefficient from positions in a time series.
%
%INPUTS
%-objs: array containing tracer information (position, IDs, ect.), for
%   example output from fo5_rp.m.
%-timestep: timestep between successive frames.
%-scale: unit conversion, if desired.
%-method: implemented options
%   -msd: Less precise. Fit line to find slope of msd curve.
%   -cve: covariance based estimator (Vestergaard et al. 2014)
%       More precise.
%-geometry: geometry of the surface on which the object is diffusing;
%       either 'planar' for a plane, or 'spherical'. Testing in progress on
%       'spherical' option.
%-rotating: also calculate Dr (rotational diffusion coefficient). Requires
%   orientation data of tracers.
%-verificationIDs: calculate chi^2 values for
%   a null hypothesis that tracers are diffusing.
%-plotopt: Display histogram of periodogram values for verification
%   purposes. The theoretical line should approximate the values of the
%   experimental histogram.
%-Rad: Radius of sphere for spherical geometry option
%-mindt: minimum timestep to consider (msd method only).
%-cutoff: Length of tracks to consider. See msdanalyze_rp.m. (msd method
%   only)
%
% OUTPUTS
%- Dt: Translational diffusion coefficients, organized as
%       Dt(1,:) = estimate
%       Dt(2,:) = uncertainty in estimate
%       Dt(3,:) = number of frames in track
%       Dt(4,:) = particle localization error
%       Dt(5,:) = Signal-to-Noise ratio
%       Dt(6,:) = track ID
%   and columns are different tracks.
%- Dr: rotational diffusion coefficient, same formatting
%- chi2: chi^2 value for a null hypothesis of pure diffusion
%
% Tristan Hormel

%% Pre-define some variables with commonly used values

% Note that you don't need to input variables that irrelevant for your
% chosen method.

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = 0;
end

if ~exist('mindt', 'var') || isempty(mindt)
    mindt = 1;
end

if ~exist('cutoff', 'var') || isempty(mindt)
    cutoff = .25;
end

if ~exist('verificationIDs', 'var') || isempty(verificationIDs)
    verificationIDs = [];
end

if ~exist('geometry', 'var') || isempty(geometry)
    geometry = 'planar';
end

if ~exist('Rad', 'var') || isempty(Rad)
    Rad = NaN;
end


% R is a constant related to the camera shutter. If the shutter is open for
% the entire exposure time, R = 1/6. It is hard to see why we wouldn't want
% this to be the case, so I'm hard coding the value. If you close the
% shutter for part of the exposure, you will need to modify this.

R = 1/6;

% some sanity checks

if le(timestep, 0)
    error('timestep should be positive and greater than 0'); 
end

if le(scale, 0)
    error('scale should be positive and greater than 0');
end

if ~strcmp(method,'msd') && ~strcmp(method,'cve')
    error('method should be cve or msd'); 
end

if ~strcmp(geometry,'spherical') && ~strcmp(geometry,'planar')
    error('method should be planar or spherical');  
end

if ~strcmp(geometry,'spherical') && le(Rad,0)
    error('Radius should be positive when using a spherical geometry');
end

%% calculations and function calls

% See functions at bottom

switch lower(method)
    
    case 'msd'
        
        [Dr, Dt] = msdEstimate(msad, msd, cutoff, plotopt, mindt, rotating);
        
    case 'cve'
        
        [Dr, Dt] = cveEstimate(objs, scale, timestep, rotating, geometry, Rad, R);
        
end

chi2 = verifyDiffusion(objs,verificationIDs,scale,plotopt,geometry,timestep,R);

end

%% cve estimate

function [Dr, Dt] = cveEstimate(objs, scale, timestep, rotating, geometry, Rad, R)

utrk = unique(objs(6,:));  % all the unique track ids

% pre-allocate variables.

D1 = zeros(1,numel(utrk)); D2 = zeros(1,numel(utrk)); sig1 = zeros(1,numel(utrk));
sig2 = zeros(1,numel(utrk)); epsilon1 = zeros(1,numel(utrk));
epsilon2 = zeros(1,numel(utrk)); sigD1 = zeros(1,numel(utrk));
sigD2 = zeros(1,numel(utrk)); Dt = zeros(6,numel(utrk)); Dr = zeros(6,numel(utrk));


for j = 1:numel(utrk)
    
    trtmp = objs(:, ismember(objs(6,:),utrk(j))); % data for jth tracer
    
    N = size(trtmp,2)-1;
    
    switch lower(geometry)
        
        case 'planar'
            
            delta1 = diff(trtmp(1,:))*scale; % displacements in x
            delta2 = diff(trtmp(2,:))*scale; % displacements in y
            
        case 'spherical'
            
%             delta1 = 2*Rad*scale.*asin(sqrt(sin(diff(trtmp(1,:))./2).^2)); % polar angular displacements
%             delta2 = 2*Rad*scale.*asin(sqrt(cos(trtmp(1,1:end-1)).*cos(trtmp(1,2:end)).*sin(diff(trtmp(2,:))./2).^2)); % azimuthal angular displacements
%             
            delta1 = diff(trtmp(1,:))*Rad*scale; % displacements in x
            delta2 = diff(trtmp(2,:))*Rad*scale; % displacements in y
% 
%             delta1 = scale*Rad*diff(trtmp(1,:));
%             delta2 = scale*Rad*sin(trtmp(1,1:end-1)).*diff(trtmp(2,:));
% 
    end
    
    % componentwise cve estimate. Equations 14 - 18 in Vestergaard et al.
    
    D1(j) = mean(delta1.^2)/2/timestep + mean(delta1(1:(end-1)).*delta1(2:end))/timestep; % estimate for 1st coord.
    D2(j) = mean(delta2.^2)/2/timestep + mean(delta2(1:(end-1)).*delta2(2:end))/timestep; % estimate for 2nd coord.
    sig1(j) = sqrt(R*mean(delta1.^2) + (2*R-1)*mean(delta1(1:(end-1)).*delta1(2:end))); % localization error 1st coord.
    sig2(j) = sqrt(R*mean(delta2.^2) + (2*R-1)*mean(delta2(1:(end-1)).*delta2(2:end)));
    epsilon1(j) = sig1(j)^2/D1(j)/timestep - 2*R;
    epsilon2(j) = sig2(j)^2/D2(j)/timestep - 2*R;
    sigD1(j) = D1(j)*sqrt((6+4*epsilon1(j) + 2*epsilon1(j)^2)/N+4/N/N*(1+epsilon1(j))^2); % std of estimate 1st. coord.
    sigD2(j) = D2(j)*sqrt((6+4*epsilon2(j) + 2*epsilon2(j)^2)/N+4/N/N*(1+epsilon2(j))^2);
    
    % we can average over the estimates from each component to obtain a
    % better total estimate of D
    
    Dt(1,j) = mean([D1(j) D2(j)]);
    Dt(2,j) = mean([sigD1(j) sigD2(j)]);
    Dt(3,j) = numel(delta1); %number of frames
    Dt(4,j) = mean([sig1(j) sig2(j)]);
    Dt(5,j) = sqrt(2*Dt(1,j)*timestep)/Dt(4,j); %SNR
    Dt(6,j) = utrk(j); %particle ID
    
    if rotating
        
        deltatheta = diff(trtmp(7,:));
        
        Dr(1,j) = mean(deltatheta.^2)/2/timestep + mean(deltatheta(1:(end-1)).*deltatheta(2:end))/timestep;
        Dr(2,j) = sqrt(R*mean(deltatheta.^2)+(2*R-1)*mean(deltatheta(1:(end-1)).*deltatheta(2:end)));
        Dr(3,j) = numel(deltatheta);
        Dr(4,j) = R*mean(deltatheta.^2)+(2*R-1)*mean(deltatheta(1:(end-1)).*deltatheta(2:end));
        Dr(5,j) = sqrt(2*Dr(1,j)*timestep)/Dr(4,j); %SNR
        Dr(6,j) = utrk(j); %particle ID
        
    else
        
        Dr = NaN;
        
    end
    
end

end

%% verification

function [chi2] = verifyDiffusion(objs,verificationIDs,scale,plotopt,geometry,timestep,R)

% We compare the data's periodogram with the expected spectrum for pure
% diffusion. This is an important check, and should generally be run for
% each tracer. Note that, among other things, this can tell us if tracer -
% tracer interactions are significant!

% Keeping notation consistent with Vestergaard et al. 2012. See page 6, appendix J.

% Here we are considering the periodogram of the experimental data. 
% This periodogram will deviate from the theoretical prediction by some 
% amount. The test is to see if this deviation follows its expected form, 
% which is a gamma function.

% Outline: (1) calculate experimental periodogram, (2) calculate theoretical
% periodogram, (3) form distribution of ratios, (4) calculate expected 
% distribution of ratios, (5) get chi^2 value by comparing results from
% (3) and (4).

% Also note that this process uses machinery from the cve estimate- here,
% we are ONLY checking that tracers are diffusing, NOT that the value of
% our D estimate is any good.

if isempty(verificationIDs)
    
    chi2 = NaN;
    return
end

if ~(ismember(objs(6,:),verificationIDs))
    warning('verificationIDs do not match obj IDs, ignoring input');
    chi2 = NaN;
    return
end

chi2 = zeros(1,numel(verificationIDs));

for j = 1:length(verificationIDs)
    
    % (1) calculate experimental periodogram
    
    trtmp = objs(:, ismember(objs(6,:),verificationIDs(j))); % separating out track j
    
    N = size(trtmp,2)-1; % total number of displacements in track
    
    switch lower(geometry)
        
        case 'planar'
            
            delta1 = diff(trtmp(1,:))*scale; % displacements in x
            
        case 'spherical'
            
            delta1 = diff(trtmp(1,:))*scale*Rad; % polar angular displacements
            
    end
    
    deltaxhat = zeros(1,N);
    
    for k = 1:N
        deltaxhat(k) = timestep*sum(sin(pi*k.*(1:N)/(N+1)).*delta1); % discrete sine transform of x
    end
    
    Phatexp = (2.*deltaxhat.^2)/((N+1)*timestep); % experimental periodogram
    
    % (2) get theoretical periodogram
    
    D = mean(delta1.^2)/2/timestep + mean(delta1(1:(end-1)).*delta1(2:end))/timestep; % estimate for 1st coord. from cve
    sigloc = sqrt(R*mean(delta1.^2)+(2*R-1)*mean(delta1(1:(end-1)).*delta1(2:end))); % estimate for 2nd coord. from cve
    
    Phattheory = 2*D*timestep^2 + 2*(sigloc^2*timestep-2*D*R*timestep^2).*(1-cos((pi.*(1:N))./(N+1)));
    
    % (3) get histogram of the ratio of experimental to theoretical
    % periodogram values.
    
    epsilonhatexp = Phatexp./Phattheory;
    
    % rule of thumb: # of bins ~ sqrt(# of data points). Since we should
    % have at least 100 frames for good statistics, and probably not more
    % than 1000, choose 10 bins.
    
    binwidth = .5;
    binmax = 5;
    binnumber = binmax/binwidth;
    
    bincenters = binwidth/2:binwidth:(binmax+binwidth/2); 
    binedges = eps:binwidth:(binmax+eps); %binwidth for calculating expected form- note that the gamma function diverges at x = 0 so I've added a small offset
    
    epsilondist = hist(epsilonhatexp, bincenters); % experimental distribution of values
    
    % (4) Calculate theoretical distribution.
    
    shapeparam = 1/2; % true for pure diffusion
    scaleparam = 2; % ditto
    
    fun = @(x) gampdf(x,shapeparam,scaleparam)*N;
    
    epsilonhattheory = zeros(1,numel(bincenters)-1); % theoretical distribution
    
    for k = 1:numel(binedges)-1
        epsilonhattheory(k) = integral(fun,binedges(k),binedges(k+1)); %theoretical distribution of values
    end
    
    % (5) get chi^2. This procedure has 3 constraints, so if 10
    % bins were used, to get a reduced chi^2 divide the chi2 output by 7
    
    chi2(j) = sum(((epsilondist(1:binnumber)-epsilonhattheory).^2)./epsilonhattheory);
    
    if plotopt
        
        linepoints = .01:.01:5.01; % finer spacing so that plotted line will look nice
        
        pdfline = gampdf(linepoints,shapeparam,scaleparam)*N*binwidth;
        
        figure; bar(bincenters(1:binnumber),epsilondist(1:binnumber));
        hold on
        plot(linepoints,pdfline,'m','Linewidth',2);
        axis([0 binmax+1 0 max(epsilondist+10)]);
        ylabel 'count'
        xlabel '\epsilon'
    end
    
end

end

