%% Reading S2 matrix (S2 generated in PolSARPro)
%%
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

f1 = fopen([path 's11.bin'],'rb');
f2 = fopen([path 's12.bin'],'rb');
f3 = fopen([path 's21.bin'],'rb');
f4 = fopen([path 's22.bin'],'rb');

s11_1 = fread(f1,[2*ncols nrows],'float32');
% ENVI data type =6 ==> Complex format
% Complex (2x32 bits)data currently not supported');
% Importing as double-precision instead;==>format= 'double';
s12_1 = fread(f2,[2*ncols nrows],'float32');
s21_1 = fread(f3,[2*ncols nrows],'float32');
s22_1 = fread(f4,[2*ncols nrows],'float32');


%% Note: In general, the first row will be real part and the next row will... 
%    be imaginary part in the complex data file i.e. Real and Imaginary parts will 
%    be alternating in each row. This format may differ based on the method of storing the file.
%% https://www.imageeprocessing.com/2015/03/how-to-read-image-file-or-complex-image.html

%%
s11_real = s11_1(1:2:ncols*2,:);
s11_imag = s11_1(2:2:ncols*2,:);
s11 = complex(s11_real,s11_imag);

s12_real = s12_1(1:2:ncols*2,:);
s12_imag = s12_1(2:2:ncols*2,:);
s12 = complex(s12_real,s12_imag);

s21_real = s21_1(1:2:ncols*2,:);
s21_imag = s21_1(2:2:ncols*2,:);
s21 = complex(s21_real,s21_imag);

s22_real = s22_1(1:2:ncols*2,:);
s22_imag = s22_1(2:2:ncols*2,:);
s22 = complex(s22_real,s22_imag);

fclose('all');
tic

%%
%Image visualization
f10 = figure('Name', 'Matrix element-S11');
set(gca,'FontSize',15)
imagesc(abs(s11)')
axis('image');
colormap('gray');
colorbar('FontSize', 15);
caxis([0 0.8]);
%% Testing with SNAP RS2>>Radiometric calibrate check i and q channel
% Snap pixel location column, row : java starts with 0,0
% MATLAB starts with 1,1
% Hence in MATLAB: img (col-1, row-1) == SNAPimg (col, row)
fprintf('S11_real (1385,1197) == %.4f \n', s11_real(1385,1197));
fprintf('S11_imag (1385,1197) == %.4f\n', s11_imag(1385,1197));
fprintf('S11 (1385,1197) == %.4f+%.4fj\n', real(s11(1385,1197)), imag(s11(1385,1197)));


%% Creating colour image
% Combine separate color channels into an RGB image.
dbl = abs(s11) -  abs(s22);     % Shh-Svv
vol = abs (s21); %Shv
surf = abs(s11) + abs(s22);
rgbImage = cat(3, dbl', vol', surf');

%% Enhance Composite with a Contrast Stretch
stretched_rgbImage1 = imadjust(rgbImage,stretchlim(rgbImage));
%Visualization of RGB image
f11 = figure;
set(f11,'name','PauliRGB (Contrast Stretch-lin)','numbertitle','off');
imshow(stretched_rgbImage1);
axis('image');
%caxis([-35 0]);
title('Pauli-RGB-MATLAB');

%% Enhance Composite with a Decorrelation Stretch
stretched_rgbImage2 = decorrstretch(rgbImage, 'Tol', 0.05);
%Visualization of RGB image
f12 = figure;
set(f12,'name','PauliRGB (Decorrelation Stretch)','numbertitle','off');
imshow(stretched_rgbImage2);
axis('image');
% caxis([0 1]);
title('Pauli-RGB-MATLAB2');

%%
%Save visualized RGB .png file in same path
figname_png12 = strcat([path,'PauliRGB(decorrstretch).png']);
print(f12,figname_png12,'-dpng')
fclose('all');
%%  end of code
