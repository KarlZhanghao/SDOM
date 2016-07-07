function [A, B, phy, datamat, fitmat] = omCal( img3D, ang)
%% sine fit for the whole image
% use bw mask to reduce calculation
% image data type must double
% fun = A * cos( 2*pi / 10 * xdata - phy) + B;
%% test module
% dirRoute = 'D:\baiduyun\lab\zhanghao\_SpodTool\data\01_pre\sts_20150420cdc12Yeast006_44nm\001_t000min00s/';
% load([dirRoute, 'ReconPara.mat']);
% % read image
% for kk = 1 : size( ReconPara.img, 2)
%     imgTmp = imread( [ReconPara.PreDir, ReconPara.img{kk}]);
%     img3D(:,:,kk) = im2double( imgTmp);
%     ang(kk) = ReconPara.ang(kk);
% end
%%
absAng = 0;
%% 
hh = size( img3D, 1);
ww = size( img3D, 2);
num = size( img3D, 3);
B = mean( img3D, 3);
B = B/max(max(B));
A = zeros( hh, ww);
phy = zeros( hh, ww);
level = graythresh( B);
mask = im2bw( B, level*1.1);
datamat = zeros( size( img3D));
fitmat = zeros( size( img3D));
%% 
for kk = 1 : hh
    for ll = 1 : ww
        if mask( kk,ll) == 0
            continue;
        end
        yvec = reshape( img3D( kk,ll,:), num, 1);
        [ATmp, BTmp, phyTmp, yfit] = omSinefit( yvec, ang);
        phy(kk,ll) = phyTmp;
        B(kk,ll) = BTmp;
        A(kk,ll) = ATmp;
        datamat( kk, ll, :) = yvec(:);
        fitmat (kk, ll,:) = yfit(:);
    end
end
phy = phy + absAng;
phy = mod( phy, 180);
