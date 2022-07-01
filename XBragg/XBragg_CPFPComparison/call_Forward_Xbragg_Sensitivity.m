%% Call forward X-Bragg for CPSAR and Full-pol 
%% Low roughness condition by beta = 10 degree
% @author: Dr. Dipankar Mandal
%%  ---------------------------------------------------------------------------------------
%   ---------------------------------------------------------------------------------------
%   Copyright (C) 2022 by Microwave Remote Sensing Lab, IITBombay http://www.mrslab.in
%  
%   This program is free software; you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the Free
%   Software Foundation; either version 3 of the License, or (at your option)
%   any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
%   more details.
%  
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, see http://www.gnu.org/licenses/
%   ---------------------------------------------------------------------------------------
%%
clear;

%% Incidence angle (SAR geometry)
inc = 35; %in degree
%% Volumetric soil moisture given real dielectric constant of soil
epr = 2:0.1:30;
eps = complex(epr,0.00);
[xx,yy] = size(epr);
for i = 1:yy
mv(i) = 0.01*(-5.3+(2.92.*epr(i))-(0.055.*epr(i).*epr(i))+(0.0004.*epr(i).*epr(i).*epr(i)));  %% Topp et al. 1980 model
end

%% Compact-pol
%% Plotting Sigma RH and RV at different width beta (0<b<90) distribution
b10 = 10; % width of distribution of beta in degree

for i = 1:yy
[sigma0RH10(i), sigma0RV10(i)] = xbragg_compactpol(b10,eps(i),inc);
end
%% Plotting with varying eps
figure
line(mv,sigma0RH10,'Color','black','LineStyle','-','LineWidth',1.2)
hold on
text(max(mv),max(sigma0RH10),num2str('RH'))
grid on;
xlabel('Volumetric soil moisture (mv)');
ylabel('\sigma^\circ (dB)');
% Range of alpha [-45 +45]
%ylim([0 45])
line(mv,sigma0RV10,'Color','magenta','LineStyle','-','LineWidth',1.2)
text(max(mv),max(sigma0RV10),num2str('RV'))

%% Full-pol
%% Plotting Sigma HH, VV and HV at different width beta (0<b<90) distribution
for i = 1:yy
[sigma0HH10(i), sigma0VV10(i), sigma0HV10(i)] = xbragg_fullpol(b10,eps(i),inc);
end
%% Plotting with varying eps
line(mv,sigma0HH10,'Color','red','LineStyle','-.','LineWidth',1.2)
text(max(mv),max(sigma0HH10),num2str('HH'))
line(mv,sigma0VV10,'Color','green','LineStyle','-.','LineWidth',1.2)
text(max(mv),max(sigma0VV10),num2str('VV'))
line(mv,sigma0HV10,'Color','blue','LineStyle','-.','LineWidth',1.2)
text(max(mv),max(sigma0HV10),num2str('HV'))




%% -------------------------------------------------------------------------
% %
% b20 = 20;
% for i = 1:yy
% [sigma0RH20(i), sigma0RV20(i)] = xbragg_compactpol(b20,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH20,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH20),num2str('\beta=20'))
% 
% %
% b40 = 40;
% for i = 1:yy
% [sigma0RH40(i), sigma0RV40(i)] = xbragg_compactpol(b40,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH40,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH40),num2str('\beta=40'))
% 
% %
% b60 = 60;
% for i = 1:yy
% [sigma0RH60(i), sigma0RV60(i)] = xbragg_compactpol(b60,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH60,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH60),num2str('\beta=60'))
% 
% %
% b70 = 70;
% for i = 1:yy
% [sigma0RH70(i), sigma0RV70(i)] = xbragg_compactpol(b70,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH70,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH70),num2str('\beta=70'))
% 
% %
% b80 = 80;
% for i = 1:yy
% [sigma0RH80(i), sigma0RV80(i)] = xbragg_compactpol(b80,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH80,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH80),num2str('\beta=80'))
% 
% %
% b90 = 90;
% for i = 1:yy
% [sigma0RH90(i), sigma0RV90(i)] = xbragg_compactpol(b90,eps(i),inc);
% end
% % Plotting with varying eps
% line(mv,sigma0RH90,'Color','red','LineStyle','-','LineWidth',1.2)
% text(max(mv),max(sigma0RH90),num2str('\beta=90'))

