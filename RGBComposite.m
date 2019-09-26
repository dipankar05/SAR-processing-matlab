%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating RGB Composite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File   : RGBComposite.m
% Authors  : Dipankar Mandal
% Version  : 1.0
% Creation : 09/2019
%Institute: Microwave Remote Sensing Lab (MRSLab) http://www.mrslab.in
%Indian Institute of Technology Bombay, India
%Email: dipankar_mandal@iitb.ac.in
%--------------------------------------------------------------------------------------
%%
% Sample data: 
%AIRSAR Flevoland data: https://earth.esa.int/documents/653194/658149/AIRSAR_Flevoland
%Read C3 matrix in PolSARPro
%--------------------------------------------------------------------------------------
%%
%%
%Start code
%set path of C3 folder
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

%% Read three channels
f1 = fopen([path 'FDDdouble.bin'],'rb');
f2 = fopen([path 'FDDsurface.bin'],'rb');
f3 = fopen([path 'FDDvolume.bin'],'rb');

FDDdbl = fread(f1,[ncols nrows],'float32') + ep;
FDDsurf = fread(f2,[ncols nrows],'float32') + ep;
FDDvol = fread(f3,[ncols nrows],'float32') + ep;

fclose('all');

%% Creating colour image
% Combine separate color channels into an RGB image.
rgbImage = cat(3, FDDdbl', FDDvol', FDDsurf');

%% Enhance Composite with a Contrast Stretch
stretched_rgbImage1 = imadjust(rgbImage,stretchlim(rgbImage));
%Visualization of RGB image
f11 = figure;
set(f11,'name','FDD_Powers_RGB (Contrast Stretch-lin)','numbertitle','off');
imshow(stretched_rgbImage1);
axis('image');
%caxis([-35 0]);
title('FDD Powers');

%% Enhance Composite with a Decorrelation Stretch
stretched_rgbImage2 = decorrstretch(rgbImage, 'Tol', 0.03);
%Visualization of RGB image
f12 = figure;
set(f12,'name','FDD_Powers_RGB (Decorrelation Stretch)','numbertitle','off');
imshow(stretched_rgbImage2);
axis('image');
%caxis([-35 0]);
title('FDD Powers');

%%
%Save visualized RGB .png file in same path
figname_png11 = strcat([path,'FDDDecompositionRGB(linearcontrast).png']);
print(f11,figname_png11,'-dpng')

figname_png12 = strcat([path,'FDDDecompositionRGB(decorrstretch).png']);
print(f12,figname_png12,'-dpng')
fclose('all');
%%  end of code