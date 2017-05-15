function res = acf(data_array,dim_list,pad_flag)
%
% ; NAME:
% ;               acf (auto correlation function)
% ; PURPOSE:
% ;               Computes autocorrelatino function of 1D, 2D, or 3D array
% ;               along any combination of dimensions.  Normalizes all data_array
% ;               so that the ACF equals 1 for perfect correlation, -1 for
% ;               perfect anti-correlation, and 0 for no correlation.
% ;
% ; CATEGORY:
% ;               Image Processing
% ; CALLING SEQUENCE:
% ;               res = acf( data_array, dim_list, pad_flag )
% ; INPUTS:
% ;               data_array:   The N-dimensional input data_array array.
% ;
% ;               dim_list:     a vector specifying the dimensions along
% ;                             which to operate.  Should be input as [1 2]
% ;                             or [1 3], or 3, or [1 2 3] etc.
% ;
% ;               pad_flag:     0 for no padding, 1 for padding.
% ; OUTPUTS:
% ;               res:    autocorrelation function.
% ; PROCEDURE:
% ;               simple implementation of the Fourier autocorrelation
% ;               theroem.  When padding, data_array is normalized by
% ;               autocorrelation of a tophat function, which compensates
% ;               for zeros in padded regions.
% ; NOTES:
% ; MODIFICATION HISTORY:
% ;               Written by Tommy Angelini, The University of Florida, 1/2014.
% ;

% ;
% ;       This code 'acf.m' is copyright 2014, Thomas E. Angelini.  It
% ;       should be considered 'freeware'- and may be distributed freely in
% ;       its original form when properly attributed.
%
%

%

%..............
% Squeeze out any zeros from dim_list and sort. Check to be sure that input
% data is 3D or lower.
%..............
dim_list=squeeze(dim_list);
dim_list=sort(dim_list);

if  length(dim_list) > 3
    disp('Input array has too many dimensions.  Max D is 3.')
    return
end

%   Find number of dimesions of data_array and number of dimensions along
%   which to compute ACF.

if isvector(data_array)==1
    data_array_D = 1;
else
    data_array_D=max(size(size(data_array)));
end

corr_D=max(size(dim_list));

if  corr_D > data_array_D
    disp('Cannot compute ACF along more dimesensions than data_array array has.')
    return
end

%...............................
% The code accepts data with even or odd numbers of elements along each
% dimension, but drops one element if the number is even.
%...............................

even_flag=(floor([size(data_array)]/2)./([size(data_array)]/2)==1);

if data_array_D==1
    data_array=data_array(1:end-even_flag(1));
elseif data_array_D==2
    data_array=data_array(1:end-even_flag(1),:);
    data_array=data_array(:,1:end-even_flag(2));
else
    data_array=data_array(1:end-even_flag(1),:,:);
    data_array=data_array(:,1:end-even_flag(2),:);
    data_array=data_array(:,:,1:end-even_flag(3));
end

% There are six possible values of dim_flag: 1,2,3,4,6,9.  Handle each
% separately

dim_flag = data_array_D*corr_D;

% dim_flag==1 correspond to a 1D data array.  "real" is taken because...

if dim_flag==1
    
    data_array=data_array-nanmean(data_array(:));
    data_array=data_array/nanstd(data_array(:),1);
    
    if pad_flag==1
        
        norm_array=ones(size(data_array));
        norm_array=padarray(norm_array,(length(norm_array)-1)/2);
        data_array=padarray(data_array,(length(data_array)-1)/2);
        
        corr_fun = ifft(fft(data_array).*conj(fft(data_array)));
        corr_norm = ifft(fft(norm_array).*conj(fft(norm_array)));
        
        res = real(fftshift(corr_fun./corr_norm));
        
    else
        
        res = fftshift(ifft(fft(data_array).*conj(fft(data_array))))/(length(data_array));
        
    end
    
end

% dim_flag==2 corresponds to a 2D data array to be autocorrelated along one
% dimension.

if dim_flag==2
    
    if dim_list==1
        data_array=data_array-repmat(nanmean(data_array,1),[size(data_array,1) 1]);
        data_array=data_array./repmat(nanstd(data_array,1,1),[size(data_array,1) 1]);
    else
        data_array=data_array-repmat(nanmean(data_array,2),[1 size(data_array,2)]);
        data_array=data_array./repmat(nanstd(data_array,1,2),[1 size(data_array,2)]);
    end
    
    if pad_flag==1
        
        norm_array=ones([size(data_array)]);
        if dim_list==1
            norm_array=padarray(norm_array,[(size(data_array,1)-1)/2 0]);
            data_array=padarray(data_array,[(size(data_array,1)-1)/2 0]);
        else
            norm_array=padarray(norm_array,[0 (size(data_array,2)-1)/2]);
            data_array=padarray(data_array,[0 (size(data_array,2)-1)/2]);
        end
        
        corr_fun = ifft(fft(data_array,[],dim_list)...
            .*conj(fft(data_array,[],dim_list)),[],dim_list);
        corr_norm = ifft(fft(norm_array,[],dim_list)...
            .*conj(fft(norm_array,[],dim_list)),[],dim_list);
        
        res = real(fftshift(corr_fun./corr_norm,dim_list));
        
    else
        
        res = real(fftshift(ifft(fft(data_array,[],dim_list)...
            .*conj(fft(data_array,[],dim_list)),[],dim_list),dim_list))...
            /size(data_array,dim_list);
        
    end
    
end

% dim_flag==3 corresponds to a 3D data array to be autocorrelated along one
% dimension.

if dim_flag==3
    
    if dim_list==1
        data_array=data_array-repmat(nanmean(data_array,1),[size(data_array,1) 1 1]);
        data_array=data_array./repmat(nanstd(data_array,1,1),[size(data_array,1) 1 1]);
    elseif dim_list==2
        data_array=data_array-repmat(nanmean(data_array,2),[1 size(data_array,2) 1]);
        data_array=data_array./repmat(nanstd(data_array,1,2),[1 size(data_array,2) 1]);
    else
        data_array=data_array-repmat(nanmean(data_array,3),[1 1 size(data_array,3)]);
        data_array=data_array./repmat(nanstd(data_array,1,3),[1 1 size(data_array,3)]);
    end
    
    if pad_flag==1
        
        norm_array=ones([size(data_array)]);
        if dim_list==1
            norm_array=padarray(norm_array,[(size(data_array,1)-1)/2 0 0]);
            data_array=padarray(data_array,[(size(data_array,1)-1)/2 0 0]);
        elseif dim_list==2
            norm_array=padarray(norm_array,[0 (size(data_array,2)-1)/2 0]);
            data_array=padarray(data_array,[0 (size(data_array,2)-1)/2 0]);
        else
            norm_array=padarray(norm_array,[0 0 (size(data_array,3)-1)/2]);
            data_array=padarray(data_array,[0 0 (size(data_array,3)-1)/2]);
        end
        
        corr_fun = ifft(fft(data_array,[],dim_list)...
            .*conj(fft(data_array,[],dim_list)),[],dim_list);
        corr_norm = ifft(fft(norm_array,[],dim_list)...
            .*conj(fft(norm_array,[],dim_list)),[],dim_list);
        
        res = real(fftshift(corr_fun./corr_norm,dim_list));
        
    else
        
        res = real(fftshift(ifft(fft(data_array,[],dim_list)...
            .*conj(fft(data_array,[],dim_list)),[],dim_list),dim_list))...
            /size(data_array,dim_list);
        
    end
    
end

% dim_flag==6 corresponds to a 3D data array to be autocorrelated along two
% dimensions.

if dim_flag==6
    
    if sum(dim_list)==3
        
        mean_array=repmat(nanmean(nanmean(data_array,1),2),[size(data_array,1) size(data_array,2) 1]);
        data_array=data_array-mean_array;
        std_array=repmat(sqrt(sum(sum(data_array.^2,1),2)),...
            [size(data_array,1) size(data_array,2) 1])/size(data_array,1)/size(data_array,2);
        data_array=data_array./std_array;
        
    elseif sum(dim_list)==4
        
        mean_array=repmat(nanmean(nanmean(data_array,1),3),[size(data_array,1) 1 size(data_array,3)]);
        data_array=data_array-mean_array;
        std_array=repmat(sqrt(sum(sum(data_array.^2,1),3)),...
            [size(data_array,1) 1 size(data_array,3)])/size(data_array,1)/size(data_array,3);
        data_array=data_array./std_array;
        
    else
        
        mean_array=repmat(nanmean(nanmean(data_array,2),3),[1 size(data_array,2) size(data_array,3)]);
        data_array=data_array-mean_array;
        std_array=repmat(sqrt(sum(sum(data_array.^2,2),3)),...
            [1 size(data_array,2) size(data_array,3)])/size(data_array,2)/size(data_array,3);
        data_array=data_array./std_array;
    end
    
    if pad_flag==1
        
        norm_array=ones([size(data_array)]);
        if sum(dim_list)==3
            norm_array=padarray(norm_array,[(size(data_array,1)-1)/2 (size(data_array,2)-1)/2 0]);
            data_array=padarray(data_array,[(size(data_array,1)-1)/2 (size(data_array,2)-1)/2 0]);
            
            A=fft(fft(data_array,[],1),[],2);
            B=conj(fft(fft(data_array,[],1),[],2));
            corr_fun=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],2),1),2));
            
            A=fft(fft(norm_array,[],1),[],2);
            B=conj(fft(fft(norm_array,[],1),[],2));
            corr_norm=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],2),1),2));
            
            res = corr_fun./corr_norm/((size(data_array,1)+1)/2)/((size(data_array,2)+1)/2);
            
        elseif sum(dim_list)==4
            norm_array=padarray(norm_array,[(size(data_array,1)-1)/2 0 (size(data_array,3)-1)/2]);
            data_array=padarray(data_array,[(size(data_array,1)-1)/2 0 (size(data_array,3)-1)/2]);
            
            A=fft(fft(data_array,[],1),[],3);
            B=conj(fft(fft(data_array,[],1),[],3));
            corr_fun=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],3),1),3));
            
            A=fft(fft(norm_array,[],1),[],3);
            B=conj(fft(fft(norm_array,[],1),[],3));
            corr_norm=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],3),1),3));
            
            res = corr_fun./corr_norm/((size(data_array,1)+1)/2)/((size(data_array,3)+1)/2);
            
        else
            norm_array=padarray(norm_array,[0 (size(data_array,2)-1)/2 (size(data_array,3)-1)/2]);
            data_array=padarray(data_array,[0 (size(data_array,2)-1)/2 (size(data_array,3)-1)/2]);
            
            A=fft(fft(data_array,[],2),[],3);
            B=conj(fft(fft(data_array,[],2),[],3));
            corr_fun=real(fftshift(fftshift(ifft(ifft(A.*B,[],2),[],3),2),3));
            
            A=fft(fft(norm_array,[],2),[],3);
            B=conj(fft(fft(norm_array,[],2),[],3));
            corr_norm=real(fftshift(fftshift(ifft(ifft(A.*B,[],2),[],3),2),3));
            
            res = corr_fun./corr_norm/((size(data_array,2)+1)/2)/((size(data_array,3)+1)/2);
            
        end
        
    else
        
        if sum(dim_list)==3    
            A=fft(fft(data_array,[],1),[],2);
            B=conj(fft(fft(data_array,[],1),[],2));
            res=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],2),1),2))/size(data_array,1)^2/size(data_array,2)^2;
        elseif sum(dim_list)==4
            A=fft(fft(data_array,[],1),[],3);
            B=conj(fft(fft(data_array,[],1),[],3));
            res=real(fftshift(fftshift(ifft(ifft(A.*B,[],1),[],3),1),3))/size(data_array,1)^2/size(data_array,3)^2;
        else
            A=fft(fft(data_array,[],2),[],3);
            B=conj(fft(fft(data_array,[],2),[],3));
            res=real(fftshift(fftshift(ifft(ifft(A.*B,[],2),[],3),2),3))/size(data_array,2)^2/size(data_array,3)^2;
        end
        
    end
    
end

% dim_flag==9 corresponds to a 3D data array to be autocorrelated along all
% three dimensions. dim_flag==4 corresponds to a 2D array autocorrelated
% along all two dimensions.  They can both be handled with one set of code.

if (dim_flag==9)||(dim_flag==4)
    
    data_array=data_array-nanmean(data_array(:));
    data_array=data_array/nanstd(data_array(:));
    
    if pad_flag==1
        
        norm_array=ones([size(data_array)]);
        norm_array=padarray(norm_array,([size(data_array)]-1)/2);
        data_array=padarray(data_array,([size(data_array)]-1)/2);
        
        corr_fun = ifftn(fftn(data_array).*conj(fftn(data_array)));
        corr_norm = ifftn(fftn(norm_array).*conj(fftn(norm_array)));
        
        res = real(fftshift(corr_fun./corr_norm));
        
    else
        
        res = real(fftshift(ifftn(fftn(data_array).*conj(fftn(data_array))))...
            /size(data_array,1)/size(data_array,2)/size(data_array,3));
        
    end
    
end

end