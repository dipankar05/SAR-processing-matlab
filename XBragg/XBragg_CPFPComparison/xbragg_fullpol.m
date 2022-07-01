%% X-bragg full-pol model
% Ref: Hajnsek, I., 2001. Inversion of surface parameters using
% polarimetric SAR (Doctoral dissertation, Jena, Univ., Diss., 2001).pp.167
% https://d-nb.info/963834800/34
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
function [sigma0HH, sigma0VV, sigma0HV] = xbragg_fullpol(b,epr,inc)
%% Orientation angle -distribution
%% Considering uniform distribution in 2Beta
%b = 35; %degree

%% Relative dielectric constant
% epr = 2:0.1:40; %% real part
eps = complex(epr,0.00);

%% Incidence angle (SAR geometry)
% inc = 35; %in degree
incr = deg2rad(inc);  % in radian



%% Bragg coefficients for horizontal (Rh) and vertical (Rv) polarizations
Rh = (cos(incr) - sqrt(eps - (sin(incr).*sin(incr))))/(cos(incr) + sqrt(eps - (sin(incr).*sin(incr))));
Rv = ((eps.*cos(incr)) - sqrt(eps - (sin(incr).*sin(incr))))/((eps.*cos(incr)) + sqrt(eps - (sin(incr).*sin(incr))));


%% X-Bragg coefficients
C1 = abs((Rh + Rv)).^2;    %% |Rh + Rv|^2;
C2 = (Rh + Rv).*(conj(Rh) - conj(Rv));
C3 = 0.5 * (abs((Rh - Rv)).^2);    %% 0.5*|Rh - Rv|^2;


%% T matrix of X-bragg
t1 = C1;
t2 = C2.*(sin(2.*b.*pi/180)/(2.*b.*pi/180));
t3 = 0;
t4 = C3.*(1 + (sin(4.*b.*pi/180)/(4.*b.*pi/180)));
t5 = 0;
t6 = C3.*(1 - (sin(4.*b.*pi/180)/(4.*b.*pi/180)));

%% C matrix conversion
% Co-variance matrix C2 elements
c11 = 0.5.*(t1 + t4 + conj(t2)+ t2);
c22 = 0.5.*2.*t6;
c33 = 0.5.*(t1 + t4 - conj(t2)- t2);

sigma0HH = 10.*log10(c11);
sigma0VV = 10.*log10(c33);
sigma0HV = 10.*log10(0.5*c22);

end
