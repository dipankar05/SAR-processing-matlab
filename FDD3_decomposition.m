%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Freeman-Durden 3 component decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File   : FDD3_decomposition.m
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

%span
span = c11+c22+c33;
%%


%% Intitialization
Podd = zeros(ncols,nrows);
Pdbl = zeros(ncols,nrows);
Pvol = zeros(ncols,nrows);
%% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

for ii=startj:stopj
    for jj=starti:stopi
        
        c11s = mean2(c11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22s = mean2(c22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c33s = mean2(c33(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c13s = mean2(c13(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        spanmax = max(max(span(ii-inci:ii+inci,jj-incj:jj+incj)));

        %%
        fv=3*c22s./2; %%Volume fraction
        c11s=c11s-fv;
        c33s=c33s-fv;
        c13r=real(c13s)-(fv./3);
        c13im=imag(c13s);
        
        %% Condition--Volume scattering > total power
        if ((c11s<=ep) || (c33s<= ep))
            fv=(3/8).*(c11s+c22s+c33s+ 2.*fv);
            fd=0;
            fs=0;
        else
            %% Condition--Odd bounce or surface 
            if (c13r>=0)
                alpha=-1;
                fd=(c11s.*c33s - c13r.*c13r - c13im.*c13im)./(c11s+c33s+ 2.*c13r);
                fs=c33s-fd;
                beta=sqrt((fd+c13r).*(fd+c13r) + (c13im.*c13im))./fs;
            end 
            %% Condition--Even bounce or double
            if (c13r<0)
                beta=1;
                fs=(c11s.*c33s - c13r.*c13r - c13im.*c13im)./(c11s + c33s - 2.*c13r);
                fd=c33s-fs;
                alpha = sqrt((fs-c13r).*(fs-c13r) + (c13im.*c13im))./fd;
            end
        
        %% Scattering power calculation
        Podd(ii,jj) = fs.*(1+ beta.*beta);
        Pdbl(ii,jj) = fd.*(1+ alpha.*alpha);
        Pvol(ii,jj) = (8/3).*fv;
        
        %% ensuring non-negetive powers
        if (Podd(ii,jj)< 0)
            Podd(ii,jj) = 0;
        end
        
        if (Pdbl(ii,jj)< 0)
            Pdbl(ii,jj) = 0;
        end
        
        if (Pvol(ii,jj)< 0)
            Pvol(ii,jj) = 0;
        end
        
        %% ensuring powers not to be > span
        if (Podd(ii,jj)>spanmax)
            Podd(ii,jj) = spanmax;
        end
        
        if (Pdbl(ii,jj)>spanmax)
            Pdbl(ii,jj) = spanmax;
        end
        
        if (Pvol(ii,jj)>spanmax)
            Pvol(ii,jj) = spanmax;
        end
            
       end      
    end
   fprintf('Column: %d \n',ii);
end
    
%% Visualization of powers in dB scale
f11 = figure;
set(f11,'name','Odd_dB','numbertitle','off');
imagesc(10.*log10(Podd'));
axis('image');
caxis([-35 0]);
title('Odd-bounce power');
colormap(jet);
colorbar

%% Save visualized image as .png file in same path
% figname_png = strcat([path,'FDDsurface_dB.png']);
% print(f11,figname_png,'-dpng')


fl2 = figure;
set(fl2,'name','Even_dB','numbertitle','off');
imagesc(10.*log10(Pdbl'));
axis('image');
caxis([-35 0]);
title('Even-bounce power');
colormap(jet);
colorbar

fl3 = figure;
set(fl3,'name','Volume_dB','numbertitle','off');
imagesc(10.*log10(Pvol'));
axis('image');
caxis([-35 0]);
title('Volume power');
colormap(jet);
colorbar

%%
%File Saving in same path
%Save powers in linear power scale, not in dB scale.
f_name_100 = strcat(['FDDsurface','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,Podd, 'float32');

f_name_101 = strcat(['FDDdouble','.bin']);
fileandpath_101=strcat([path, f_name_101]);
fid_101 = fopen(fileandpath_101,'wb');
fwrite(fid_101,Pdbl, 'float32');

f_name_102 = strcat(['FDDvolume','.bin']);
fileandpath_102=strcat([path, f_name_102]);
fid_102 = fopen(fileandpath_102,'wb');
fwrite(fid_102,Pvol, 'float32');

fclose('all');

%%