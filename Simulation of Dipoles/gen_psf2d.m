function psf = gen_psf2d( fwhm, pxl_size)
%% test module
% fwhm = 250;
% pxl_size = 50;
%%
% gen grid
xx = pxl_size/2:pxl_size:fwhm*4;
if mod(length(xx),2) == 0
    xx(end+1) = xx(end) + pxl_size;
end
yy = xx;
[yy, xx] = meshgrid(yy, xx); 
%
sigma = fwhm/2.3548;
x0 = xx((length(xx)+1)/2);
y0 = x0;
%
psf = exp(- ((xx-x0).^2+(yy-y0).^2) / (2*sigma^2));
psf = psf/sum(psf(:));
