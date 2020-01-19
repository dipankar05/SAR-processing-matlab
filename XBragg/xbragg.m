function [alphaC,HC] = xbragg(b)
%% Orientation angle -distribution
%% Considering uniform distribution in 2Beta
%b = 35; %degree
    
%% Dielectric constant
epr = 0.00001:1:40;
eps = complex(epr,0.0);
[x,y] = size(epr);
% AA = 0.0*ones(x,y);

for i = 1:y
    
    %%
    % % epr = 30;
    % % eps = complex(epr,0.0);
    
    %% Incidence angle (SAR geometry)
    inc = 35; %in degree
    incr = 35*pi/180;  % in radian
    
    %% Fresnel coefficients for horizontal (Rh) and vertical (Rv) polarizations
    Rh = (cos(incr) - sqrt(eps(i) - (sin(incr)*sin(incr))))./(cos(incr) + sqrt(eps(i) - (sin(incr)*sin(incr))));
    Rv = ((eps(i)*cos(incr)) - sqrt(eps(i) - (sin(incr)*sin(incr))))./((eps(i)*cos(incr)) + sqrt(eps(i) - (sin(incr)*sin(incr))));
    
    
    %% X-Bragg coefficients
    C1 = abs(Rh - Rv)^2;    %% |Rh + Rv|^2;
    C2 = (Rh + Rv)*(conj(Rh) - conj(Rv));
    C3 = abs(Rh + Rv)^2;    %% |Rh - Rv|^2;
    
    
    %% T3 elements of X-bragg
    t11 = 0.5*C1;
    t12 = 0.5*C2*(sin(2*b*pi/180)/(2*b*pi/180));
    t13 = 0;
    t21 = conj(t12);
    t22 = 0.25*C3*(1+(sin(4*b*pi/180)/(4*b*pi/180)));
    t23 = 0;
    t31 = 0;
    t32 = 0;
    t33 = 0.25*C3*(1-(sin(4*b*pi/180)/(4*b*pi/180)));
    
    T3 = [t11 t12 t13;t21 t22 t23;t31 t32 t33];
    
    %% H-Alpha decomposition
    [evec_v, eval] = eig(T3);
    
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
    
    %Lambda 2
    eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
    
    eval_norm2(eval_norm2 < 0) = 0;
    eval_norm2(eval_norm2 > 1) = 1;
    
    %Lambda 3
    eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
    
    eval_norm3(eval_norm3 < 0) = 0;
    eval_norm3(eval_norm3 > 1) = 1;
    
    
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
    alphaC(i) = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
        eval_norm3*alpha3*180./pi);
    
    %Entropy
    HC(i) = -eval_norm1*log10(eval_norm1)./log10(3) - ...
        eval_norm2*log10(eval_norm2)./log10(3) - ...
        eval_norm3*log10(eval_norm3)./log10(3);
    
end

end


