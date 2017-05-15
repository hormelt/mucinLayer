function res = gaussfit2dTH(data)

rows=1:size(data,1);
cols=1:size(data,2);

[cols rows]=meshgrid(cols,rows);

xd=cols(:); % x= first column
yd=rows(:); % y= second column
dataL=data(:);   % your data f(x,y) (in column vector)
dataW=dataL;
dataW(dataW==0)=inf;
dataW(dataW==inf)=0.99*min(dataW);

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0 0 0 0 0 0],...
    'Upper',[inf size(data,2) size(data,2) size(data,1) size(data,1) 2^12],...
    'Startpoint',[max(data(:)) size(data,2)/2 size(data,2)/2 size(data,2)/2 size(data,2)/2 0],...
    'MaxFunEvals',1E12,...
    'MaxIter',1E12,...
    'Weights',(1./dataW).^0);
f = fittype('A1*exp(-((x1-A2).^2)/(2*A3^2)).*exp(-((x2-A4).^2)/(2*A5^2))+A6','independent',{'x1','x2'},'options',s);

[c2,gof2]=fit([xd yd],dataL,f,'Exclude',dataL < min(dataW));

% fitted_surface = c2.A1*exp(-((cols-c2.A2).^2)/(2*c2.A3^2)).*exp(-((rows-c2.A4).^2)/(2*c2.A3^2))+c2.A5;

errors = confint(c2);
errorbars = (errors(2,:)-errors(1,:))/2;

res = [c2.A1 errorbars(1) c2.A2 errorbars(2) c2.A3 errorbars(3) c2.A4 errorbars(4) c2.A5 errorbars(5) c2.A6 errorbars(6)];

% figure;
% imagesc(data);
% hold on

% Y=AT.*x.^BT;
fitted_surface = c2.A1*exp(-((cols-c2.A2).^2)/(2*c2.A3^2)).*exp(-((rows-c2.A4).^2)/(2*c2.A5^2))+c2.A6;
% surface(x1,x2,fitted_surface);

% imagesc([data fitted_surface 2*(data-fitted_surface)]);

% hold off

end