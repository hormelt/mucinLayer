clear

frmend=22;
zend = 19;

zaxis = linspace(-6.75,6.75,19);

rng('shuffle');

x = round((rand(1000,1)*199+1));
y = round((rand(1000,1)*199+1));

stack=zeros(512,512,zend,frmend);

h = zeros(length(x),frmend);
herror = h;

A0out = h;
A1out = h;
d2out = h;
z2out = h;

A0error = h;
A1error = h;
d2error = h;
z2error = h;

A1T = 800;
A0T = 100;
d1T = 0.5;
d2T = 0.5;
z1T = -6;
z2T = 1;

for frame = 1:frmend
    frame
    for z = 1:zend
        
        
        temp=double(imread(['image_tristan_001_t',num2str(frame,'%03g'),'_z',num2str(z,'%03g'),'_c001.tif']));
        %         stack(:,:,z,frame) = temp((256-100):(256+99),(256-100):(256+99));
        stack(:,:,z,frame) = temp;
        
    end
    
end


for step = 1:length(x)
    step
    xstep = x(step);
    ystep = y(step);
    
    for frame = 1:frmend
        
        A1T = 800;
        A0T = 100;
        d2T = 0.5;
        z2T = 1;
        
        data = squeeze(stack(ystep,xstep,:,frame));
        
        subdata = squeeze(stack(ystep,xstep,5:end,frame));
        subz = zaxis(5:end);
        
        % check guesses by plotting
%                 I = A1T*(1-erf(d2T*(zaxis-z2T)))+A0T;
%         
%                                 hold off
%                                 plot(zaxis,data)
%                                 hold on
%         
%                                 plot(zaxis,I,'r')
        
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0 0 0 -2],...
            'Upper',[1e3 2e3 1000 5],...
            'Startpoint',[A0T A1T d2T z2T],...
            'MaxFunEvals',1E8,...
            'MaxIter',1E8);%,...
        %                 'Weights',(max(subdata(:))./subdata).^0.5);
        f = fittype('A1*(1-erf(d2*(x-z2)))+A0','independent','x','options',s); %1-erf(d2*(x-z2)) is probability that d2*(x-z2) is greater than sqrt(2)*sigma*d2*(x-z2) from mean. 4 parameter fit to 1000 data points
        [c2,gof2]=fit(subz',subdata,f);
        c2;
        
        A0T=c2.A0;
        A1T=c2.A1;
        d2T=c2.d2;
        z2T=c2.z2;
        
        I = A1T*(1-erf(d2T*(zaxis-z2T)))+A0T;
%         
%         hold off
%         plot(zaxis,data)
%         hold on
%         
%         plot(zaxis,I,'r')
%         getframe
        
        h(step,frame)=z2T-zaxis(2);
        A0out(step,frame)=A0T;
        A1out(step,frame)=A1T;
        d2out(step,frame)=d2T;
        z2out(step,frame)=z2T;
        
        errors=confint(c2);
        
        herror(step,frame)=(errors(2,4)-errors(1,4))/2+0.75/2;
        A0error(step,frame)=(errors(2,1)-errors(1,1))/2;
        A1error(step,frame)=(errors(2,2)-errors(1,2))/2;
        d2error(step,frame)=(errors(2,3)-errors(1,3))/2;
        z2error(step,frame)=(errors(2,4)-errors(1,4))/2;
        
    end
end



csvwrite('hrt.csv',h);
csvwrite('hrt_error.csv',herror);

csvwrite('A0rt.csv',A0out);
csvwrite('A0rt_error.csv',A0error);
csvwrite('A1rt.csv',A1out);
csvwrite('A1rt_error.csv',A1error);
csvwrite('d2rt.csv',d2out);
csvwrite('d2rt_error.csv',d2error);
csvwrite('z2rt.csv',z2out);
csvwrite('z2rt_error.csv',z2error);

csvwrite('xsample.csv',x);
csvwrite('ysample.csv',y);
