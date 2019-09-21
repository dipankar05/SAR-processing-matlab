%%
%H-Alpha decomposition and segmentation with H-Alpha plane ftom T3 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File   : Eigen_decomposition.m
% Authors  : Dipankar Mandal
% Version  : 1.0
% Creation : 09/2019
%Institute: Microwave Remote Sensing Lab (MRSLab) http://www.mrslab.in
%Indian Institute of Technology Bombay, India
%Email: dipankar_mandal@iitb.ac.in
%--------------------------------------------------------------------------------------
%%
%% EIGENDECOMPOSITION OF COHERENCY MATRIX
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

f1 = fopen([path 'T11.bin'],'rb');
f2 = fopen([path 'T12_real.bin'],'rb');
f3 = fopen([path 'T12_imag.bin'],'rb');
f4 = fopen([path 'T13_real.bin'],'rb');
f5 = fopen([path 'T13_imag.bin'],'rb');
f6 = fopen([path 'T22.bin'],'rb');
f7 = fopen([path 'T23_real.bin'],'rb');
f8 = fopen([path 'T23_imag.bin'],'rb');
f9 = fopen([path 'T33.bin'],'rb');

t11_T1 = fread(f1,[ncols nrows],'float32') + ep;
t12_T1 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
t21_T1 = conj(t12_T1);
t13_T1 = complex( fread(f4,[ncols nrows],'float32') , fread(f5,[ncols nrows],'float32')) + ep;
t31_T1 = conj(t13_T1);
t22_T1 = fread(f6,[ncols nrows],'float32') + ep;
t23_T1 = complex( fread(f7,[ncols nrows],'float32') , fread(f8,[ncols nrows],'float32')) + ep;
t32_T1 = conj(t23_T1);
t33_T1 = fread(f9,[ncols nrows],'float32') + ep;

fclose('all');
tic

%% Intitialization

B1 = zeros(ncols,nrows);
B2 = zeros(ncols,nrows);
B3 = zeros(ncols,nrows);
B4 = zeros(ncols,nrows);
B5 = zeros(ncols,nrows);
B6 = zeros(ncols,nrows);
B7 = zeros(ncols,nrows);
B8 = zeros(ncols,nrows);
B9 = zeros(ncols,nrows);

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
        
        t11s = mean2(t11_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t12s = mean2(t12_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t13s = mean2(t13_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        t21s = mean2(t21_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t22s = mean2(t22_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t23s = mean2(t23_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        t31s = mean2(t31_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t32s = mean2(t32_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t33s = mean2(t33_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        T_T1 = [t11s t12s t13s; t21s t22s t23s; t31s t32s t33s];
        
        [evec_v, eval] = eig(T_T1);
        
        %% Eigenvalues
        
        eval_diag = (sort(diag(eval)))';
        
        if (eval_diag(1) < 0)
            eval_diag(1) = 0;
        end
        
        if (eval_diag(2) < 0)
            eval_diag(2) = 0;
        end
        
        if (eval_diag(3) < 0)
            eval_diag(3) = 0;
        end
        
        %Lambda 1
        eval_norm1 = (eval_diag(3))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm1(eval_norm1 < 0) = 0;
        eval_norm1(eval_norm1 > 1) = 1;
        
        B1(ii,jj) = eval_norm1;
        
        %Lambda 2
        eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
        B2(ii,jj) = eval_norm2;
        
        %Lambda 3
        eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm3(eval_norm3 < 0) = 0;
        eval_norm3(eval_norm3 > 1) = 1;
        
        B3(ii,jj) = eval_norm3;
        
        %% Eigenvectors
        
        %Alpha 1
        eig_vec_r1 = real(evec_v(1,3));
        eig_vec_c1 = imag(evec_v(1,3));
        
        alpha1 = acos(sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1));
        B4(ii,jj) = alpha1*180./pi;
        
        %Alpha 2
        eig_vec_r2 = real(evec_v(1,2));
        eig_vec_c2 = imag(evec_v(1,2));
        
        alpha2 = acos(sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2));
        B5(ii,jj) = alpha2*180./pi;
        
        %Alpha 3
        eig_vec_r3 = real(evec_v(1,1));
        eig_vec_c3 = imag(evec_v(1,1));
        
        alpha3 = acos(sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3));
        B6(ii,jj) = alpha3*180./pi;
        
        %Cloude Alpha
        B7(ii,jj) = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
            eval_norm3*alpha3*180./pi);
        
        %Entropy
        B8(ii,jj) = -eval_norm1*log10(eval_norm1)./log10(3) - ...
            eval_norm2*log10(eval_norm2)./log10(3) - ...
            eval_norm3*log10(eval_norm3)./log10(3);
        
        %Anisotropy
        B9(ii,jj) = (eval_norm2 - eval_norm3)./(eval_norm2 + eval_norm3);
        
    end
    fprintf('Column: %d \n',ii);
end

%% Plot of H/Alpha without zones and segmentation
f1 = figure;
set(f1,'name','Alpha/Entropy','numbertitle','off');
plot(B8(:), B7(:),'.', 'MarkerSize',0.2);
grid on;
xlabel('Entropy');
ylabel('Alpha');

%% Vusualisation du resultat (pas obligatoire)
fl1 = figure;
set(fl1,'name','Lambda 1','numbertitle','off');
imagesc(B1');
axis('image');
caxis([0 1]);
title('Lambda 1');
colormap(jet);
colorbar

fl2 = figure;
set(fl2,'name','Lambda 2','numbertitle','off');
imagesc(B2');
axis('image');
caxis([0 1]);
title('Lambda 2');
colormap(jet);
colorbar

fl3 = figure;
set(fl2,'name','Lambda 3','numbertitle','off');
imagesc(B3');
axis('image');
caxis([0 1]);
title('Lambda 3');
colormap(jet);
colorbar

f4 = figure;
set(f4,'name','Alpha 1','numbertitle','off');
imagesc(B4');
axis('image');
caxis([0 90]);
colormap(jet);
title('Alpha 1');
colorbar

f5 = figure;
set(f5,'name','Alpha 2','numbertitle','off');
imagesc(B5');
axis('image');
caxis([0 90]);
colormap(jet);
title('Alpha 2');
colorbar

f6 = figure;
set(f6,'name','Alpha 3','numbertitle','off');
imagesc(B6');
axis('image');
caxis([0 90]);
colormap(jet);
title('Alpha 3');
colorbar

f7 = figure;
set(f7,'name','Cloude Alpha','numbertitle','off');
imagesc(B7');
axis('image');
caxis([0 90]);
colormap(jet);
title('Cloude Alpha');
colorbar

f8 = figure;
set(f8,'name','Entropy','numbertitle','off');
imagesc(B8');
axis('image');
caxis([0 1]);
colormap(jet);
title('Entropy');
colorbar

f9 = figure;
set(f9,'name','Anisotropy','numbertitle','off');
imagesc(B9');
axis('image');
caxis([0 1]);
colormap(jet);
title('Anisotropy');
colorbar
%% 

%% H-Alpha segmentation
thre_mat = ones(ncols,nrows);
ent = B8;
theta = B7;
for i = 1:ncols
   for j = 1:nrows
       if ent(i,j)>=0 && ent(i,j) <= 0.5
           if theta(i,j) >= 0 && theta(i,j) <42.5
               thre_mat(i,j) = 9;
           end
           if theta(i,j) >=42.5 && theta(i,j)<47.5
               thre_mat(i,j) = 8;
           end
           if theta(i,j) >=47.5 && theta(i,j) <=90
               thre_mat(i,j) = 7;
           end
       end
       % % % % % % % % 
       if ent(i,j)>0.5 && ent(i,j) <=0.9
            if theta(i,j) >= 0 && theta(i,j) <40
               thre_mat(i,j) = 6;
           end
           if theta(i,j) >=40 && theta(i,j)<50
               thre_mat(i,j) = 5;
           end
           if theta(i,j) >=50 && theta(i,j) <=90
               thre_mat(i,j) = 4;
           end
       end
       % % % % % % % % 
       if ent(i,j)>0.9 && ent(i,j) <=1.0
           if theta(i,j) >= 0 && theta(i,j) <40
               thre_mat(i,j) = 3;
           end
           if theta(i,j) >=40 && theta(i,j)<60
               thre_mat(i,j) = 2;
           end
           if theta(i,j) >=60 && theta(i,j) <=90
               thre_mat(i,j) = 1;
           end
       end
     % % % % % % % %
   end
   fprintf('H-Alpha segmetation started...Column: %d \n',i);
end
clusHAlpha = thre_mat;



%%
% plotting H-Alpha clusters in plane
X = [B8(:) B7(:)];  % [H Alpha]
classID = clusHAlpha;
classNames = {'Z1','Z2', 'Z3', 'Z4', 'Z5','Z6', 'Z7', 'Z8', 'Z9'}; %# one name per class
%Creating colour pallete
colors = [125 125 125; 255 150 150; 105 255 105; 255 80 80; 80 255 80; 80 80 255; 154 0 0; 0 154 0; 0 0 154]./255;

%% Plotting final figure
%%
f100 = figure('Name', 'H-Alpha');
% set(gca,'FontSize',15)
hold on
for i=1:9 %# there are 9 classes
id = classID == i;
plot(X(id,1),X(id,2),'.','Color',colors(i,:),'MarkerSize',10,'DisplayName',classNames{i})
end
%legend('show')

%% Zone boundaries
axis([0 1 0 90])
% line([x1 x2],[y1 y2])
hold on
plot([0.5 0.5],[0 90])
hold on
line([0.9 0.9],[0 90])
hold on
line([0 0.5],[42.5 42.5])
hold on
line([0 0.5],[47.5 47.5])
hold on
line([0.5 0.9],[40 40])
hold on
line([0.5 0.9],[50 50])
hold on
line([0.9 1.0],[40 40])
hold on
line([0.9 1.0],[60 60])
%-------------------------------------------------------------------------
% zone tagging
text(0.95,75,{'Z1'});
text(0.95,50,{'Z2'});
text(0.95,20,{'Z3'});
text(0.70,75,{'Z4'});
text(0.70,45,{'Z5'});
text(0.70,20,{'Z6'});
text(0.25,75,{'Z7'});
text(0.25,44,{'Z8'});
text(0.25,20,{'Z9'});
%-----------------------------------
grid on;
xlabel('Entropy');
ylabel('Alpha');
%%  
%%
%Save visualized H-Alpha plane as .png file in same path
figname_png = strcat([path,'HAlphaplane.png']);
print(f100,figname_png,'-dpng')
fclose('all');
%%  end of code

