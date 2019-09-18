%%
%Generating Stokes parameters for 2x2 RH-RV OR LH-LV Covariance matrix
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
%% Read C2 matrix
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

%% C2 matrix
f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C22.bin'],'rb');

c11 = fread(f1,[ncols nrows],'float32') + ep;
c12 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c21 = conj(c12);
c22 = fread(f4,[ncols nrows],'float32') + ep;

fclose('all');


%% Window processing
%with optional window operation: keep window=1 for pixel based

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing
%%

%Simulate circular polarization
chi_transmit = input('Input Tau (Ellipticity: +ve for LC and -ve RC transmit || range -45 to +45 degree): '); % RC = -theta and LC = +theta


%% Intitialization
cpr = zeros(ncols,nrows);
dop = zeros(ncols,nrows);
%s1 = zeros(ncols,nrows);%for saving s1 or s2..

for ii=startj:stopj
    for jj=starti:stopi
        
        %% C2 matrix
        
        c11c = mean2(c11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12c = mean2(c12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c21c = mean2(c21(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22c = mean2(c22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        C0 = [c11c c12c; c21c c22c];
        
        %% LH-LV/RH-RV
        
        % Stokes Parameter
        
        s0 = c11c + c22c;
        s1 = c11c - c22c;
        %s1(ii,jj) = c11c - c22c; %if want to save the file
        s2 = c12c + c21c;
        
        if (chi_transmit >= 0)
            s3 = (1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
        end
        
        if (chi_transmit < 0)
            s3 = -(1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
        end
        
        %% Stokes vector child products
        
        SC = (s0 - s3)./2;
        OC = (s0 + s3)./2;
        
        %%Degree of polarization
        dop(ii,jj) = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        
        %% CPR
        cpr(ii,jj) = SC./OC;

        %%
        span = trace(C0);

        %% Eigen values
        e_v = sort(eig(C0),'descend');
        e_v1 = e_v(1);%first eigen value
        e_v2 = e_v(2);%Second eigen value


    end
    fprintf('Column: %d \n',ii);
end

%% Visualization
imagesc(dop')