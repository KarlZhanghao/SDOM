clear all; clc; warning off;
%% Type parameters
%
% parameters for reconstruction of 500 nm beads
% DataSet = 'Bead500nm';
% Ang = [42;26;10;174;158;142;126;110;94;78;62];
% ReconPara = [0.06; 10; 2001];
% adjR2th = 0.6;
% vecZoom = 2;
%
% parameters for reconstruction of Septin data
% DataSet = 'Septin1';
% Ang = [155;139;123;107;92;76;60;44;28;12;176;160;144;128;112;97;81;65;49;33;17;1];
% ReconPara = [0.075; 20; 201];
% adjR2th = 0.4;
% vecZoom = 2;
%
% parameters for reconstruction of Septin data
% DataSet = 'Septin2';
% Ang = [79;68;57;47;36;25;15;4;173;163;152];
% ReconPara = [0.075; 10; 201];
% adjR2th = 0.6;
% vecZoom = 3;
%
% parameters for reconstruction of Neuronal Spine data
DataSet = 'Spine1';
Ang = [0;18;36;54;72;90;108;126;144;162];
ReconPara = [0.05;10;351];
adjR2th = 0.4;
vecZoom = 3;
%
% parameters for reconstruction of Neuronal Spine data
% DataSet = 'Spine2';
% Ang = [0;18;36;54;72;90;108;126;144;162];
% ReconPara = [0.075;10;251];
% adjR2th = 0.4;
% vecZoom = 3;
%% Change Working Directory and Read Data
disp( ['Start Processing Dataset: ', DataSet])
disp( 'Reading Data ...')
% Change Working Directory
s = what;
curDir = [s.path, '/'];
cd( curDir);
addpath('./Util/');
% Read and Display Image Data 
dataDir = strrep(fullfile(curDir, '../Data/', DataSet, '/'), '\', '/');
wfDir = [dataDir, 'WideField/'];
psfDir = [dataDir, 'PSF/'];
info = dir( wfDir);
for kk = 3 : length( info)
    img(:,:,kk-2) = imread( [wfDir, info(kk).name]);
end
img = double(img);
imgShow = sum( img, 3);
figure(1)
imshow( imgShow, [], 'InitialMagnification','fit');
title('Wide Field fluorescent image');
% Read and Display PSF data
info = dir( psfDir);
psf = imread( [psfDir, info(3).name]);
psf = double(psf);
psf = psf / sum(psf(:));
figure(2)
imshow( psf, [], 'InitialMagnification','fit');
title( 'Point Spread Funtion of the System');
%% Orientation Mapping of conventional wide field images
disp( 'Orientation Mapping of conventional wide field images ...')
% A cropped image of small area is recommended for orientation mapping
% imgFit = img(RectArea(3):RectArea(4), RectArea(1):RectArea(2), :);
imgFit = img;
% Create directory to store OM-WideField data
omWfDir = [dataDir, 'OM_WideField/'];
mkdir( omWfDir)
% OM-WideField calculation
[A, B, phy, datamat, fitmat] = omCal( imgFit, Ang);
% calculate normalized RMSE
for kk = 1 : size( datamat, 3)
    Arep(:,:,kk) = A;
    Brep(:,:,kk) = B;
end
normFit = (fitmat - Brep)./Arep; normFit( Arep==0) = 0;
normData = (datamat - Brep)./Arep; normData( Arep==0) = 0;
% normRMSE = sqrt( sum((normData-normFit).^2,3) / size( datamat,3));
normTV = sum( normData.^2, 3);
normEV = sum((normData-normFit).^2,3);
adjR2 = 1-normEV./normTV; adjR2(normTV==0) = 0;
mask = adjR2>adjR2th;
% mask = normRMSE<RMSEth;
% Display results
[x,y] = meshgrid( 1:size(datamat,2), 1:size(datamat, 1));
figure(3)
hold off
imshow( sum(imgFit,3), [], 'InitialMagnification','fit')
colormap('Hot')
hold on
OUF = A ./ B; OUF(OUF>1) = 1; OUF(mask==0) = 0;
maxOUF = max(OUF(:));
v1 = OUF.*cos(phy/180*pi); v1(mask==0) = 0; 
u1 = OUF.*(sin(phy/180*pi)); u1(mask==0) = 0;
quiver(x,y,v1,u1,vecZoom*maxOUF, 'color', 'b', 'LineStyle', '-');
v2 = OUF.*cos((phy+180)/180*pi); v2(mask==0) = 0; 
u2 = OUF.*(sin((phy+180)/180*pi)); u2(mask==0) = 0;
quiver(x,y,v2,u2,vecZoom*maxOUF, 'color', 'b', 'LineStyle', '-');
% save data
print( 3, '-dtiff', '-r800', [omWfDir, 'OM-WF.tif']);
imwrite( OUF, [omWfDir, 'OUF-WF.tif']);
save([omWfDir, 'OM-WF_check.mat'], 'A', 'B', 'Ang', 'datamat', 'fitmat', 'phy', 'OUF', 'u1', 'u2', 'v1', 'v2', 'adjR2');
%% super resolution reconstrution using SDOM
% save data and parameters to run python reconstruction program
tmpDir = [dataDir, 'tmp/'];
mkdir( tmpDir)
data.ReconPara = ReconPara;
data.image = imgFit;
data.psf = psf;
save( [tmpDir, 'ReconData.mat'], 'data');
% run reconstrution
disp( 'Reconstructing super-resolved images ...')
omSdomDir = [dataDir, 'OM_SDOM/'];
mkdir( omSdomDir)
command = ['python ', 'SDOM.py ', dataDir];
[status, cmdout] = system( command);
%
load( [tmpDir, 'ReconResult.mat']);
rmdir( tmpDir, 's');
img = data;
imgSave = sum( img, 3);
imgSave = uint16( imgSave / max(imgSave(:)) * 65535) ;
imwrite( imgSave, [omSdomDir, 'SR_SDOM.tif'])
%% Orientation Mapping of SDOM images
disp( 'Orientation Mapping of SDOM images ...')
% OM-WideField calculation
[A, B, phy, datamat, fitmat] = omCal( img, Ang);
% calculate normalized RMSE
for kk = 1 : size( datamat, 3)
    Arep(:,:,kk) = A;
    Brep(:,:,kk) = B;
end
normFit = (fitmat - Brep)./Arep; normFit( Arep==0) = 0;
normData = (datamat - Brep)./Arep; normData( Arep==0) = 0;
normTV = sum( normData.^2, 3);
normEV = sum((normData-normFit).^2,3);
adjR2 = 1-normEV./normTV; adjR2(normTV==0) = 0;
mask = adjR2>adjR2th;
% Display results
[x,y] = meshgrid( 1:size(datamat,2), 1:size(datamat, 1));
figure(4)
hold off
imshow( sum(img,3), [], 'InitialMagnification','fit')
colormap('Hot')
hold on
OUF = A ./ B; OUF(OUF>1) = 1; OUF(mask==0) = 0;
maxOUF = max(OUF(:));
v1 = OUF.*cos(phy/180*pi); v1(mask==0) = 0; 
u1 = OUF.*(sin(phy/180*pi)); u1(mask==0) = 0;
quiver(x,y,v1,u1,maxOUF, 'color', 'b', 'LineStyle', '-');
v2 = OUF.*cos((phy+180)/180*pi); v2(mask==0) = 0; 
u2 = OUF.*(sin((phy+180)/180*pi)); u2(mask==0) = 0;
quiver(x,y,v2,u2,maxOUF, 'color', 'b', 'LineStyle', '-');
% save data
print( 4, '-dtiff', '-r800', [omSdomDir, 'OM-SDOM.tif']);
imwrite( OUF, [omSdomDir, 'OUF-SDOM.tif']);
save([omSdomDir, 'OM-SDOM_check.mat'], 'A', 'B', 'Ang', 'datamat', 'fitmat', 'phy', 'OUF', 'u1', 'u2', 'v1', 'v2', 'adjR2');
