
%% Transmitted Stokes
%%V polarization
chitr = 0;
psitr = 90;

% %%H polarization
% chitr = 0;
% psitr = 0;

% %%Right circular -R
% chitr = -45;
% psitr = 0;

% %%Left circular -L
% chitr = 45;
% psitr = 0;

Str = [1; cos(2*psitr*pi/180)*cos(2*chitr*pi/180); sin(2*psitr*pi/180)*cos(2*chitr*pi/180); sin(2*chitr*pi/180)];

%% Kennaugh Matrix of Elementary targets        
K_d = [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]; %dihedral
K_nd = [0.625 0.375 0 0; 0.375 0.625 0 0; 0 0 -0.5 0; 0 0 0 0.5];  %Narrow dihedral
K_t = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]; %trihedral
K_c = [0.625 0.375 0 0; 0.375 0.625 0 0; 0 0 0.5 0; 0 0 0 -0.5]; %Cylinder
K_lh = [1 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1]; %left helix
K_rh = [1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]; %right helix
K_rv = [1 0 0 0; 0 0.5 0 0; 0 0 0.5 0; 0 0 0 0]; %45 deg rotated dihedral
K_id = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; %Ideal depolarizer
%% Received Stokes
Src = K_t*Str; %% Change K_x as per elementary target

%%Properties of received wave
m = sqrt(Src(2)^2 + Src(3)^2 + Src(4)^2)/Src(1); %degree of polarization
chir = 0.5 * (180/pi) * asin((1)*Src(4))/(m*Src(1)) % Ellipticity

%%Normalized Stokes
Srcn = Src/Src(1);

%% Dummy Stokes from C2
% Srcc = [0.0166; 0.0065; 0.0038; -0.0018];
% Srcn = Srcc/Srcc(1);

%% Create and plot a sphere with a radius equal to 1. 
C = [0 0 0]; % center vector (h,k,l)
r = 1; % radius
n = 50; % n^2 number of faces
[x,y,z] = sphere(n); % creates x, y, and z matrices
x = r*x + C(1); % adjusts x matrix with r and C info
y = r*y + C(2); % adjusts y matrix with r and C info
z = r*z + C(3); % adjusts z matrix with r and C info
a = round((n+1)/2); % creates an index from n info

axes1 = axes('Parent',figure);
% s = surf(x,y,z,'facecolor','blue','facealpha',0.10,'edgecolor','blue',...
%     'edgealpha',0.5,'linewidth',0.05); % surface plot of the sphere
s = surf(x,y,z,'facealpha',0.35); % surface plot of the sphere
s.EdgeColor = 'none';
shading interp 
light               % create a light
lighting gouraud    % preferred method for lighting curved surfaces
material dull    % set material to be dull, no specular highlights
axis equal
hold on
plot3(x(a,:),y(a,:),z(a,:),'b','linewidth',1.5); % plot of great circle
plot3([C(1),Srcn(2)],[C(2),Srcn(3)],[C(3),Srcn(4)],'k','linewidth',1.0); % plot of segment between C and A
plot3([0,0],[0,0],[0,1],'k','linewidth',1.0);
plot3([0,1],[0,0],[0,0],'k','linewidth',1.0);
plot3([0,0],[0,1],[0,0],'k','linewidth',1.0);
xlabel('S_2')
ylabel('S_3')
zlabel('S_4')
%% Plot elementary target
plot3(Srcn(2),Srcn(3),Srcn(4),'or','markersize',7,'markerfacecolor','red'); % plot of point A
view(axes1,[125 20]);
