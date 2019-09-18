%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading C3 matrix of quad-pol SAR data (preprocessed in PolSARPro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File   : Reading_C3matrix.m
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
%Generate C3 or T3 matrix in PolSARPro
%--------------------------------------------------------------------------------------
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

%% C3 matrix
f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C13_real.bin'],'rb');
f5 = fopen([path 'C13_imag.bin'],'rb');
f6 = fopen([path 'C22.bin'],'rb');
f7 = fopen([path 'C23_real.bin'],'rb');
f8 = fopen([path 'C23_imag.bin'],'rb');
f9 = fopen([path 'C33.bin'],'rb');

c11 = fread(f1,[ncols nrows],'float32') + ep;
c12 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c13 = complex( fread(f4,[ncols nrows],'float32') , fread(f5,[ncols nrows],'float32')) + ep;
c21 = conj(c12);
c22 = fread(f6,[ncols nrows],'float32') + ep;
c23 = complex( fread(f7,[ncols nrows],'float32') , fread(f8,[ncols nrows],'float32')) + ep;
c31 = conj(c13);
c32 = conj(c23);
c33 = fread(f9,[ncols nrows],'float32') + ep;

fclose('all');

%%
%Image visualization
f10 = figure('Name', 'Matrix element-C11');
set(gca,'FontSize',15)
imagesc(c11')
axis('image');
colormap('gray');
colorbar('FontSize', 15);
caxis([0 0.05]);

%%
%Save visualized image as .png file in same path
figname_png = strcat([path,'C11.png']);
print(f10,figname_png,'-dpng')

%%
%File Saving in same path
f_name_100 = strcat(['C11_mod','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,c11, 'float32');

fclose('all');

%%