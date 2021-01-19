
chi_in = 0;
c11 = 0.011554;
c12_real = 0.001905;
c12_imag = 0.000878;
c12 = complex(c12_real,c12_imag);
c21 = conj(c12);
c22 = 0.005065;

%% LH-LV/RH-RV
% Stokes Parameter
s0 = c11 + c22;
s1 = c11 - c22;
s2 = c12 + c21;
        
if (chi_in >= 0)
    s3 = (1i.*(c12 - c21)); %The sign is according to RC or LC sign !!
end
if (chi_in < 0)
    s3 = -(1i.*(c12 - c21)); %The sign is according to RC or LC sign !!
end

dop = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0)
DOCP = (s3)./(dop.*s0); 
Chir = 0.5.*((180/pi).*asin(DOCP))
        
%% m-chi decomposition     
sin_2chi = (s3./(dop.*s0));     
surface_m_chi = ((dop.*(s0).*(1+sin_2chi)./2))
double_bounce_m_chi = ((dop.*(s0).*(1-sin_2chi)./2))
diffused_m_chi = (((s0).*(1-dop)))
