
%% H-Alpha curve
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




%% Plotting with different beta 0<b<90
%%
b1 = 10;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=10'))

%%
b1 = 20;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=20'))

%%
b1 = 40;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=40'))

%%
b1 = 60;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=60'))

%%
b1 = 70;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=70'))

%%
b1 = 80;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=80'))

%%
b1 = 90;
[alphaC1,HC1] = xbragg(b1);
%% Plotting with varying eps
line(HC1,alphaC1,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC1),max(alphaC1),num2str('\beta=90'))
