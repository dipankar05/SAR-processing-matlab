clear 
%% need HA_Anisotropytheta function
m = 0:0.01:1;
[x,y] = size(m);
%% Curve -1 --Lower curve
%log3(9) = log10(9) / log10(3)

for i = 1:y
H1(i) = ((-1)/(1+2*m(i)))*(log10(power(m(i),2*m(i))/power((1+2*m(i)),(1+2*m(i))))/log10(3));
a1(i) = (m(i)*180)/(2*m(i)+1);
end
figure
line(H1,a1)
hold on
%% Curve -II
%%limit-1 0<m<0.5
m = 0:0.01:0.5;
[x,y] = size(m); 
for i = 1:y
H21(i) = ((-1)/(1+2*m(i)))*(log10(power(2*m(i),2*m(i))/power((1+2*m(i)),(1+2*m(i))))/log10(3));
a21(i) = 180/2;
end
line(H21,a21) 
 
%%limit-1 0.5<m<1
m = 0.5:0.01:1;
[x,y] = size(m); 
for i = 1:y
H22(i) = ((-1)/(1+2*m(i)))*(log10(power((2*m(i)-1),(2*m(i)-1))/power((1+2*m(i)),(1+2*m(i))))/log10 (3));
a22(i) = 180/(1+2*m(i));
end
line(H22,a22) 

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
%%
%%
%%
%%
%% Plotting the anisotrpic media loci with fixed A, varying theta
%% A = 0.0
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.0*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0'))
%% A = 0.20
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.20*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0.2'))
%% A = 0.40
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.40*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0.4'))
%% A = 0.60
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.60*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0.6'))
%% A = 0.80
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.80*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0.8'))
%% A = 0.90
theta = 0.00001:1:90;
[x,y] = size(theta);
AA = 0.90*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('A=0.9'))

%%
%%
%%
%%
%% Plotting the anisotrpic media loci with fixed theta, varying A
%% theta = 10 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 10*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
%%https://in.mathworks.com/help/matlab/ref/line.html
text(max(HC),max(alphaC),num2str('T=10'))
% 
%% theta = 20 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 20*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
%%https://in.mathworks.com/help/matlab/ref/line.html
text(max(HC),max(alphaC),num2str('T=20'))

%% theta = 30 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 30*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('T=30'))

%% theta = 50 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 50*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('T=50'))

%% theta = 70 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 70*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('T=70'))

%% theta = 90 deg
AA = 0.0001:0.01:1;
[x,y] = size(AA);
theta = 90*ones(x,y);

for i = 1:y
[HC(i),alphaC(i)] = HA_Anisotropytheta(AA(i),theta(i));
end
line(HC,alphaC,'Color','blue','LineStyle','-.','LineWidth',1.2)
text(max(HC),max(alphaC),num2str('T=90'))
