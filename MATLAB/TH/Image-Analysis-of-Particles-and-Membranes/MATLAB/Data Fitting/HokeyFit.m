function [ etafit, afit, chifit ] = HokeyFit( Dr, Dt, T )
% Calculate viscosities from rotational and translational diffusion coefficients.
%
% INPUTS
%   -Dr : array of rotational diffusion coeffiecients from msadtr.m
%   -Dt : array of translational diffusion coefficients from msdtr_rp.m
%   -T : temperature, in Kelvin
%
% OUTPUTS
%   -etafit : best fit for 2D viscosity
%   -chifit : chi^2 value of best fit
% Tristan Hormel
% last modified 11/2/12 - include fits using logarithms

%constant definitions
kB = 1.38e-23;  % Boltzmann's constant, J/K
etaw = 9e-4;  % viscosity of water, Pa.s
gamma = 0.577; % Euler's constant
T = kB*T;
alpha = 100e-9;

%initial ranges of viscosity to check
eta = logspace(-10,-5,1000);

Dtval = Dt(1,:)*1e-12;  % translational diffusion coefficients, m^2/s
sigDt = Dt(2,:)*1e-12; % uncertainty in Dt, m^2/s
Drval = Dr(1,:);  % rotational diffusion coefficients, radians^2/s
sigDr = Dr(2,:);  % uncertainty in Dr
nanind = find(isnan(Dtval) | isnan(sigDt) | isnan(Drval) | isnan(sigDr));
Dtval(nanind) = [];
sigDt(nanind) = [];
Drval(nanind) = [];
sigDr(nanind) = [];

%formatting to avoid looping
eta = repmat(eta,length(Drval),1)';
Dtval = repmat(Dtval,length(eta(:,1)),1);
Drval = repmat(Drval,length(eta(:,1)),1);
sigDt = repmat(sigDt,length(eta(:,1)),1);
sigDr = repmat(sigDr,length(eta(:,1)),1);


%%
%Chi^2 minimization, several methods

chi2rough = sum((Dtval - DtHokey(etaw, alpha, Drval, eta, T)).^2./sigDt.^2,2);
etafitindrough = find(chi2rough==min(chi2rough)); %find the minimum chi^2 value
etafitrough = eta(etafitindrough,1);
%Next we construct an array containing measured values of Dr and values to
%their immediate left and right. We will use these points to
%calculate the slope of the HPW contour at each Dr value. This
%slope will then be able to convert Dr error into Dt error.
Drvalarray = [];
for j = 1:size(Drval,2)
    Drvalarray = [Drvalarray Drval(1,j)-.01*Drval(1,j) Drval(1,j) Drval(1,j)+.01*Drval(1,j)];
end
Dtvalarray = DtHokey(etaw, alpha, Drval, eta, T);
%We now have Dr and Dt values from which to determine the slope at
%the measured Dr values using HPW, so let's do it. Sort of wierd
%indexing here- all that is happening is that we are calculating
%the slope using the neighboring points.
for j = 1:size(Drval,2)
    neighborinds = (3*j-2):(3*j);
    [inter(j) slope(j)] = fitline(Drvalarray(neighborinds), Dtvalarray(neighborinds));
end
%From this slope we can convert Dr error to Dt error. The next step
%is just formatting.
slope = repmat(slope,length(eta(:,1)),1);
w = 1./(sigDt.^2+sigDr.^2.*slope.^2);

chi2 = sum((Dtval - DtHokey(etaw, alpha, Drval, eta, T)).^2.*w,2);
        
etafitind = find(chi2==min(chi2)); %find the minimum chi^2 value
etafit = eta(etafitind,1);
chifit = chi2(etafitind);

chi2uncert = abs(chi2-2*chifit);
lowboundind = find(chi2uncert==min(chi2uncert(1:etafitind)));
highboundind = find(chi2uncert==min(chi2uncert(etafitind:end)));
etalowbound = eta(lowboundind,1);
etahighbound = eta(highboundind,1);

etafit = [etafit, abs(etahighbound-etalowbound)];
afit = 1;

%% Figures

disp('temporarily turning off negative numbers on log. axis warning.')
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');


%Plot contours to see validity of fit visually
disp('temporarily turning off negative numbers on log axis warning.')
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');

%Plot contours to see validity of fit

figure('name', 'Diffusion Measurements and Fit');
errorxy([Dr(1,:)' Dt(1,:)' Dr(2,:)' Dt(2,:)'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
    'EdgeColor', [0.0 0.3 0.7], 'FaceColor', [0.3 0.6 0.9], 'MarkSize',8,'WidthEB',1,'ScaleX', 'log');
ax = axis;
Drmin = ax(1); Drmax = ax(2); %Draw the contour across the entire plot
Drpoints = logspace(log10(Drmin),log10(Drmax),100); %Use these points to calculate the contour

Fitcontour = DtHokey(etaw, alpha, Drpoints, etafit(1), T); %the contour line
semilogx(Drpoints,Fitcontour*1e12,'Color',[.4 .7 .3],'LineWidth',2);
errorxy([Dr(1,:)' Dt(1,:)' Dr(2,:)' Dt(2,:)'], 'ColX', 1, 'ColY', 2, 'ColXe', 3, 'ColYe', 4, ...
    'EdgeColor', [0.0 0.3 0.7], 'FaceColor', [0.3 0.6 0.9], 'MarkSize',8,'WidthEB',1,'ScaleX', 'log');


warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');

%% functions

    function Dtcalc = DtHokey(etaw, alpha, Dr, eta, T)
        X = log(eta)-log(etaw)-1/2.*log(T/4/pi./eta./Dr-2*etaw*alpha^3./eta);
        Dtcalc = T.*X./(2.*eta-3*etaw*alpha.*X)/2/pi;
    end

end

