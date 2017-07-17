clear variables;
%% define samples
% system PSF
fwhm = 250;
pxl_size = 50;
psf = gen_psf2d(250, 50);
% Poarization of excitations
theta = (0 : 18 : 179) / 180 * pi;
% dipoles in the sample
gs = zeros( 21, 73);
g = zeros( size(gs,1), size(gs,2), size(theta,2));
f = [11,9,0;11,11,0;11,27,0;11,29,0;11,45,0;11,47,0;11,63,0;11,65,0];
ang0 = rand(1)*180;
f(1:2:8,3) = ang0/180*pi;
f(2:8:8,3) = (ang0-30)/180*pi;
f(4:8:8,3) = (ang0-50)/180*pi;
f(6:8:8,3) = (ang0-70)/180*pi;
f(8:8:8,3) = (ang0-90)/180*pi;
%% photons emitted from the dipole under polarized excitation
for kk = 1 : size(f, 1)
    gs( f(kk,1), f(kk,2)) = 1;
end
imwrite( gs, 'sample.tif')
for kk = 1 : size( theta, 2)
    for ll = 1 : size( f, 1)
        gTmp = zeros( size(gs));
        gTmp(f(ll,1), f(ll,2)) = cos(theta(kk)-f(ll,3))^2;
        g(:,:,kk) = g(:,:,kk) + gTmp;
    end
    imwrite( g(:,:,kk), ['sample\sample_', num2str(kk), '.tif'])
end
g = g * 65535;
g_ave = sum(sum(sum(g))) / size(g,1) / size(g,2)/ size(g, 3);
lamb1 = 1 / g_ave;
%% generate Gaussian noise
b = randn( size( gs));
b = b/std(b(:));
b = 100 + 30*b;
% b_ = dct2(b);
% b_ave = sum(sum(abs(b_)))/size(b_,1)/size(b_,2);
% lamb2 = 1 / b_ave;
%% generated convolved images acquired by the camera
u = zeros( size(gs,1), size(gs,2), size(theta,2));
img = zeros( size(gs,1), size(gs,2), size(theta,2));
for kk = 1 : size( theta, 2)
    u(:,:,kk) = conv2( g(:,:,kk), psf, 'same') + b;
    img(:,:,kk) = imnoise( uint16(u(:,:,kk)), 'poisson');
    img = uint16(img);
    imwrite( img(:,:,kk), ['data\', num2str(kk-1), '.tif'])
end
