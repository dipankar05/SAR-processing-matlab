function [HC,alphaC] = HA_Anisotropytheta(A,theta)
% A = 0.2;
% theta = 20;
%% running the function
%[HC,alphaC] = HA_Anisotropytheta(0.20,25)

%%
th = theta*pi/180;

%%
t11 = (1+A)^2;
t12 = (1-A*A)*sin(2*th)/(2*th);
t13 = 0;
t21 = t12';
t22 = ((A-1)*(A-1))*((4*th)+sin(4*th))/(8*th);
t23 = 0;
t31 = 0;
t32 = 0;
t33 = ((A-1)*(A-1))*((4*th)-sin(4*th))/(8*th);


T = [t11 t12 t13; t21 t22 t23; t31 t32 t33];
[evec_v, eval] = eig(T);



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

p1 = eval_norm1;

%Lambda 2
eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));

eval_norm2(eval_norm2 < 0) = 0;
eval_norm2(eval_norm2 > 1) = 1;

p2 = eval_norm2;

%Lambda 3
eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));

eval_norm3(eval_norm3 < 0) = 0;
eval_norm3(eval_norm3 > 1) = 1;

p3 = eval_norm3;

%% Eigenvectors

%Alpha 1
eig_vec_r1 = real(evec_v(1,3));
eig_vec_c1 = imag(evec_v(1,3));

alpha1 = acos(sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1));


%Alpha 2
eig_vec_r2 = real(evec_v(1,2));
eig_vec_c2 = imag(evec_v(1,2));

alpha2 = acos(sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2));


%Alpha 3
eig_vec_r3 = real(evec_v(1,1));
eig_vec_c3 = imag(evec_v(1,1));

alpha3 = acos(sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3));


%Cloude Alpha
alphaC = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
    eval_norm3*alpha3*180./pi);

%Entropy
HC = -eval_norm1*log10(eval_norm1)./log10(3) - ...
    eval_norm2*log10(eval_norm2)./log10(3) - ...
    eval_norm3*log10(eval_norm3)./log10(3);

%Anisotropy
% AC = (eval_norm2 - eval_norm3)./(eval_norm2 + eval_norm3);

fprintf('Entropy: %f, Alpha: %f \n',HC,alphaC);
end