function [ res ] = gaussFit2D( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



rows = 1:size(data,1);
cols = 1:size(data,2);

[cols rows] = meshgrid(cols,rows);

xd=cols(:);
yd=rows(:);

% I think that probably this only makes sense for logaritmic data, which
% I'm not doing. Also modify code so that I won't use threshholded values

dataW=data; % Set weights for fitting- why are we weighting toward low values?
% dataW(dataW==0)=inf; %preventing weighting of 0 pixels from fucking everything up- probably a better way of doing this
% dataW(dataW==inf)=0.99*min(dataW(:));

s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0 0 0 0 0],...
        'Upper',[1e4 size(data,2) size(data,2) size(data,1) 100],...
        'Startpoint',[max(data(:)) size(data,2)/2 size(data,2)/2 size(data,2)/2 0],...
        'MaxFunEvals',1E12,...
        'MaxIter',1E12);
    
f = fittype('A1*exp(-((x1-A2).^2)/(2*A3^2)).*exp(-((x2-A4).^2)/(2*A3^2))+A5','independent',{'x1','x2'},'options',s);
[c2,gof2]=fit([xd yd],data(:),f,'Exclude',data(:)<min(dataW(:)));
fitted_surface = c2.A1*exp(-((cols-c2.A2).^2)/(2*c2.A3^2)).*exp(-((rows-c2.A4).^2)/(2*c2.A3^2))+c2.A5;

if sum(dataW(:)>0)>6

    errors = confint(c2);
    errorbars = (errors(2,:)-errors(1,:))/2;
    
    res = [c2.A1 errorbars(1) c2.A2 errorbars(2) c2.A3 errorbars(3) c2.A4 errorbars(4) c2.A5 errorbars(5)];
    
else
    
    res = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    
end

end

