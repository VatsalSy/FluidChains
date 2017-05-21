%% takes in the data from analytical.mat saved by analyticalVatsalpart1.m
%% then an analytical model is tested for the given case. 
clc
clear
close
load('analytical.mat')
plot(xs,-zs,'k-', 'Linewidth', 3)
hold on
plot(xs,zs,'k-', 'Linewidth', 3)
hold on
set(streamslice(X/d, Z/d, UW, WW, 1,'arrows','nearest'),'Color', 0.75*[0 0 1], 'Linewidth', 3);
hold on

%% Flow parameters
Fr = 2.5;
alpha = 60/2;  % angle in degree
d = 0.005;
g = 9.8067;

%% Fiting parameters
o = 2; %offset
eta = 10;

%% rejection parameters
xrej = -10;

%% calculation part
uj =  Fr*sqrt(9.8067*d); % velocity in m/s
x = linspace(0,max(xs)+o,100);
em = [0.25 0.5 0.75 1]; % coefficient of restitution

for j = 1:1:length(em)
e = em(j);
psi0 = atand(e*tand(alpha));
b = eta*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
yd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
y = zeros(length(x),1);
for jint = 2:1:length(x)
    y(jint) = trapz(x(1:jint)*d,yd(1:jint))/d;
end
x1 = linspace(min(xs) - 0.5 , 0 , 10);
y1 = -tand(alpha)*x1;
xp = [x1 x];
yp = [y1 y'];
xp = xp - o;
%% Rejection Part
xptemp = zeros(length(xp),1);
yptemp = zeros(length(xp),1);
counter = 0;
for count = 1:1:length(xp)
     if xp(count) < xrej
         continue;
     else
             counter = counter + 1;
             xptemp(counter) = xp(count);
             yptemp(counter) = yp(count);
     end
end
xp = xptemp(1:counter);
yp = yptemp(1:counter);

%% Final plotting
plot(xp,yp,'g.',xp,-yp,'r.','markers',12)
hold on
end
camup([-1 0 0])
x0 = linspace(min(xs),max(xs),100);
plot(x0(1:2:end),zeros(length(x)/2,1),'k.','markers',12)
plot(x0(2:2:end),zeros(length(x)/2,1),'r.','markers',12)
axis equal tight
xlabel('X/d_j')
ylabel('Y/d_j')
