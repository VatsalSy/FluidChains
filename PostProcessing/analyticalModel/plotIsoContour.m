clc
clear
close
error = [0 0 0 0];
subplot(1,5,1)
%% One
place = 'case1-18/sim0.250.gfs';
d = 0.005;
alpha = 30;
Ldomain = 10*d;
g = 9.8067;
xmin = -0.1*Ldomain;
xc = -0.077*Ldomain + 0.5*d/sind(alpha);
xmax = 1.1*Ldomain;

% Fiting parameters
o = 2; %offset
eta = 2.56;


% Flow parameters
Fr = 2.5;

uj =  Fr*sqrt(9.8067*d); % velocity in m/s


% geometrical parameters
zmin = 0;
zmax = Ldomain/8;
nx = 1200;
nz = 150;
x = linspace(xmin,xmax,nx);
y = 0;
z = linspace(zmin,zmax,nz);
tic
gridfile='cartgrid3.dat';
disp('saving the 3d grid');
[X,Z,Y] = meshgrid(x,z,y);
loc=[X(:),Y(:),Z(:)];
save(gridfile,'loc','-ASCII','-SINGLE');
toc
T = structuredData(place, gridfile, X, Z, 'T');
X = X-xc;
figure
[c,~] = contour(X,Z,T,[0.5 0.5]);
close
xs = c(1,2:c(2,1)+1)/d;
zs = c(2,2:c(2,1)+1)/d;
xrej = 0;
% check and plot
x = linspace(0, (xmax - xc)/d +o,nx);
e = 1; % coefficient of restitution
psi0 = atand(e*tand(alpha));
b = eta*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
z = zeros(length(x),1);
for jint = 2:1:length(x)
    z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
end
xp = x - o;
zp = z;
% Rejection Part
xptemp = zeros(length(xp),1);
zptemp = zeros(length(xp),1);
counter = 0;
for count = 1:1:length(xp)
     if xp(count) < xrej
         continue;
     else
         counter = counter + 1;
         xptemp(counter) = xp(count);
         zptemp(counter) = zp(count);
     end
end
xp = xptemp(1:counter);
zp = zptemp(1:counter);
% analyzing the L1 error norm
xstemp = zeros(length(xs),1);
zstemp = zeros(length(xs),1);
counter = 0;
for count = 1:1:length(xs)
     if xs(count) < xrej
         continue;
     else
         counter = counter + 1;
         xstemp(counter) = xs(count);
         zstemp(counter) = zs(count);
     end
end
xs = xstemp(1:counter);
zs = zstemp(1:counter);
hf = createFitInter(xs,zs);
zs = hf(xp);
% plot(xp,zs,'r-','LineWidth',3)
% hold on
% plot(xp,-zs,'r-','LineWidth',3)
% hold on
error(1) = sum(abs((zp-zs)./zs))/length(zs);
% plot(xp(1:40:end),zp(1:40:end),'k*','MarkerSize',15)
% hold on
% plot(xp(1:40:end),-zp(1:40:end),'k*','MarkerSize',15)
% axis equal tight
% camup([-1 0 0])

yp = 2*ones(1,length(xp));
plot3(xp(1:30:end),yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
hold on
plot3(xp(1:30:end),yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
plot3(xp(1:30:end),-yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
hold on
plot3(xp(1:30:end),-yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
hold on
load('case1-18/chain.mat');
Cgly1 = smooth3((Y.^2+Z.^2)/(16*d)^2);
X = X- xc;
pgly1 = patch(isosurface(X/d, Y/d, Z/d, T, 0.5,Cgly1));
isonormals(X/d,Y/d,Z/d,T,pgly1)
pgly1.FaceColor = 'inTerp';
pgly1.EdgeColor = 'none';
daspect([1 1 1])
axis(volumebounds(X/d,Y/d,Z/d,U,V,W))
grid on
colormap(flipud(gray))
camlight
lighting gouraud

clear T U V W X Y Z
xlabel X/d_j
zlabel Z/d_j
grid on
axis equal tight
view(0,0)
camup([-1 0 0])

%% uncomment these lines to plot more iso-surfaces! Will require other "palce" as well

% subplot(1,5,4)
% %% two
% place = 'case4-31/sim0.250.gfs';
% d = 0.005;
% alpha = 30;
% Ldomain = 10*d;
% g = 9.8067;
% xmin = -0.05*Ldomain;
% xc = -0.029*Ldomain + 0.5*d/sind(alpha);
% xmax = 0.76*Ldomain;
% 
% % Fiting parameters
% o = 2; %offset
% eta = 2.377;
% 
% 
% % Flow parameters
% Fr = 2;
% 
% uj =  Fr*sqrt(9.8067*d); % velocity in m/s
% 
% 
% % geometrical parameters
% zmin = 0;
% zmax = Ldomain/8;
% nx = 1200;
% nz = 150;
% x = linspace(xmin,xmax,nx);
% y = 0;
% z = linspace(zmin,zmax,nz);
% tic
% gridfile='cartgrid3.dat';
% disp('saving the 3d grid');
% [X,Z,Y] = meshgrid(x,z,y);
% loc=[X(:),Y(:),Z(:)];
% save(gridfile,'loc','-ASCII','-SINGLE');
% toc
% T = structuredData(place, gridfile, X, Z, 'T');
% X = X-xc;
% figure
% [c,~] = contour(X,Z,T,[0.5 0.5]);
% close
% xs = c(1,2:c(2,1)+1)/d;
% zs = c(2,2:c(2,1)+1)/d;
% xrej = 0;
% 
% % check and plot
% x = linspace(0, (xmax - xc)/d +o,nx);
% e = 1; % coefficient of restitution
% psi0 = atand(e*tand(alpha));
% b = eta*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
% a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
% zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
% z = zeros(length(x),1);
% for jint = 2:1:length(x)
%     z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
% end
% xp = x - o;
% zp = z;
% % Rejection Part
% xptemp = zeros(length(xp),1);
% zptemp = zeros(length(xp),1);
% counter = 0;
% for count = 1:1:length(xp)
%      if xp(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xptemp(counter) = xp(count);
%          zptemp(counter) = zp(count);
%      end
% end
% xp = xptemp(1:counter);
% zp = zptemp(1:counter);
% % analyzing the L1 error norm
% xstemp = zeros(length(xs),1);
% zstemp = zeros(length(xs),1);
% counter = 0;
% for count = 1:1:length(xs)
%      if xs(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xstemp(counter) = xs(count);
%          zstemp(counter) = zs(count);
%      end
% end
% xs = xstemp(1:counter);
% zs = zstemp(1:counter);
% hf = createFitInter(xs,zs);
% zs = hf(xp);
% % plot(xp,zs,'r-','LineWidth',3)
% % hold on
% % plot(xp,-zs,'r-','LineWidth',3)
% % hold on
% error(2) = sum(abs((zp-zs)./zs))/length(zs);
% % plot(xp(1:40:end),zp(1:40:end),'k*','MarkerSize',15)
% % hold on
% % plot(xp(1:40:end),-zp(1:40:end),'k*','MarkerSize',15)
% % axis equal tight
% % camup([-1 0 0])
% 
% yp = 2*ones(1,length(xp));
% plot3(xp(1:30:end),yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% plot3(xp(1:30:end),-yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),-yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% alpha = 30;
% xc = -0.029*Ldomain + 0.5*d/sind(alpha);
% load('case4-31/chain.mat');
% Cgly1 = smooth3((Y.^2+Z.^2)/(16*d)^2);
% X = X- xc;
% pgly1 = patch(isosurface(X/d, Y/d, Z/d, T, 0.5,Cgly1));
% isonormals(X/d,Y/d,Z/d,T,pgly1)
% pgly1.FaceColor = 'inTerp';
% pgly1.EdgeColor = 'none';
% daspect([1 1 1])
% axis(volumebounds(X/d,Y/d,Z/d,U,V,W))
% grid on
% colormap(flipud(gray))
% camlight
% lighting gouraud
% 
% clear T U V W X Y Z
% xlabel X/d_j
% zlabel Z/d_j
% grid on
% axis equal tight
% view(0,0)
% camup([-1 0 0])

%% 
% subplot(1,5,5)
% %% three
% place = 'case5-47/sim0.250.gfs';
% d = 0.005;
% alpha = 44.5;
% Ldomain = 10*d;
% g = 9.8067;
% xmin = -0.37*Ldomain;
% xc = -0.327*Ldomain + 0.5*d/sind(alpha);
% xmax = 0.385*Ldomain;
% 
% % Fiting parameters
% o = 2; %offset
% eta = 2.513;
% 
% 
% % Flow parameters
% Fr = 2;
% 
% uj =  Fr*sqrt(9.8067*d); % velocity in m/s
% 
% 
% % geometrical parameters
% zmin = 0;
% zmax = Ldomain/8;
% nx = 1200;
% nz = 150;
% x = linspace(xmin,xmax,nx);
% y = 0;
% z = linspace(zmin,zmax,nz);
% tic
% gridfile='cartgrid3.dat';
% disp('saving the 3d grid');
% [X,Z,Y] = meshgrid(x,z,y);
% loc=[X(:),Y(:),Z(:)];
% save(gridfile,'loc','-ASCII','-SINGLE');
% toc
% T = structuredData(place, gridfile, X, Z, 'T');
% X = X-xc;
% figure
% [c,~] = contour(X,Z,T,[0.5 0.5]);
% close
% xs = c(1,2:c(2,1)+1)/d;
% zs = c(2,2:c(2,1)+1)/d;
% xrej = 0;
% 
% % check and plot
% x = linspace(0, (xmax - xc)/d +o,nx);
% e = 1; % coefficient of restitution
% psi0 = atand(e*tand(alpha));
% b = eta*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
% a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
% zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
% z = zeros(length(x),1);
% for jint = 2:1:length(x)
%     z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
% end
% xp = x - o;
% zp = z;
% % Rejection Part
% xptemp = zeros(length(xp),1);
% zptemp = zeros(length(xp),1);
% counter = 0;
% for count = 1:1:length(xp)
%      if xp(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xptemp(counter) = xp(count);
%          zptemp(counter) = zp(count);
%      end
% end
% xp = xptemp(1:counter);
% zp = zptemp(1:counter);
% % analyzing the L1 error norm
% xstemp = zeros(length(xs),1);
% zstemp = zeros(length(xs),1);
% counter = 0;
% for count = 1:1:length(xs)
%      if xs(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xstemp(counter) = xs(count);
%          zstemp(counter) = zs(count);
%      end
% end
% xs = xstemp(1:counter);
% zs = zstemp(1:counter);
% hf = createFitInter(xs,zs);
% zs = hf(xp);
% % plot(xp,zs,'r-','LineWidth',3)
% % hold on
% % plot(xp,-zs,'r-','LineWidth',3)
% % hold on
% error(3) = sum(abs((zp-zs)./zs))/length(zp);
% % plot(xp(1:40:end),zp(1:40:end),'k*','MarkerSize',15)
% % hold on
% % plot(xp(1:40:end),-zp(1:40:end),'k*','MarkerSize',15)
% % axis equal tight
% % camup([-1 0 0])
% 
% yp = 2*ones(1,length(xp));
% plot3(xp(1:30:end),yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% plot3(xp(1:30:end),-yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),-yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% load('case5-47/chain.mat');
% Cgly1 = smooth3((Y.^2+Z.^2)/(16*d)^2);
% X = X- xc;
% pgly1 = patch(isosurface(X/d, Y/d, Z/d, T, 0.5,Cgly1));
% isonormals(X/d,Y/d,Z/d,T,pgly1)
% pgly1.FaceColor = 'inTerp';
% pgly1.EdgeColor = 'none';
% daspect([1 1 1])
% axis(volumebounds(X/d,Y/d,Z/d,U,V,W))
% grid on
% colormap(flipud(gray))
% camlight
% lighting gouraud
% 
% clear T U V W X Y Z
% xlabel X/d_j
% zlabel Z/d_j
% grid on
% axis equal tight
% view(0,0)
% camup([-1 0 0])
% 

%%
% subplot(1,5,2)
% %% four
% place = 'case6-117/sim0.15.gfs';
% d = 0.005;
% alpha = 30;
% Ldomain = 10*d;
% g = 9.8067;
% xmin = -0.1*Ldomain;
% xc = -0.077*Ldomain + 0.5*d/sind(alpha);
% xmax = 0.82*Ldomain;
% 
% % Fiting parameters
% o = 2; %offset
% eta = 2.815;
% 
% 
% % Flow parameters
% Fr = 2.5;
% 
% uj =  Fr*sqrt(9.8067*d); % velocity in m/s
% 
% 
% % geometrical parameters
% zmin = 0;
% zmax = Ldomain/8;
% nx = 1200;
% nz = 150;
% x = linspace(xmin,xmax,nx);
% y = 0;
% z = linspace(zmin,zmax,nz);
% tic
% gridfile='cartgrid3.dat';
% disp('saving the 3d grid');
% [X,Z,Y] = meshgrid(x,z,y);
% loc=[X(:),Y(:),Z(:)];
% save(gridfile,'loc','-ASCII','-SINGLE');
% toc
% T = structuredData(place, gridfile, X, Z, 'T');
% X = X-xc;
% figure
% [c,h] = contour(X,Z,T,[0.5 0.5]);
% close
% xs = c(1,2:c(2,1)+1)/d;
% zs = c(2,2:c(2,1)+1)/d;
% xrej = 0;
% 
% % check and plot
% x = linspace(0, (xmax - xc)/d +o,nx);
% e = 1; % coefficient of restitution
% psi0 = atand(e*tand(alpha));
% b = eta*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
% a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
% zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
% z = zeros(length(x),1);
% for jint = 2:1:length(x)
%     z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
% end
% xp = x - o;
% zp = z;
% % Rejection Part
% xptemp = zeros(length(xp),1);
% zptemp = zeros(length(xp),1);
% counter = 0;
% for count = 1:1:length(xp)
%      if xp(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xptemp(counter) = xp(count);
%          zptemp(counter) = zp(count);
%      end
% end
% xp = xptemp(1:counter);
% zp = zptemp(1:counter);
% % analyzing the L1 error norm
% xstemp = zeros(length(xs),1);
% zstemp = zeros(length(xs),1);
% counter = 0;
% for count = 1:1:length(xs)
%      if xs(count) < xrej
%          continue;
%      else
%          counter = counter + 1;
%          xstemp(counter) = xs(count);
%          zstemp(counter) = zs(count);
%      end
% end
% xs = xstemp(1:counter);
% zs = zstemp(1:counter);
% hf = createFitInter(xs,zs);
% zs = hf(xp);
% % plot(xp,zs,'r-','LineWidth',3)
% % hold on
% % plot(xp,-zs,'r-','LineWidth',3)
% % hold on
% error(4) = sum(abs((zp-zs)./zs))/length(zp);
% % plot(xp(1:40:end),zp(1:40:end),'k*','MarkerSize',15)
% % hold on
% % plot(xp(1:40:end),-zp(1:40:end),'k*','MarkerSize',15)
% % axis equal tight
% % camup([-1 0 0])
% 
% yp = 2*ones(1,length(xp));
% plot3(xp(1:30:end),yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% plot3(xp(1:30:end),-yp(1:30:end),zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% plot3(xp(1:30:end),-yp(1:30:end),-zp(1:30:end),'r*','LineWidth',3,'MarkerSize',10)
% hold on
% load('case6-117/chain.mat');
% Cgly1 = smooth3((Y.^2+Z.^2)/(16*d)^2);
% X = X- xc;
% pgly1 = patch(isosurface(X/d, Y/d, Z/d, T, 0.5,Cgly1));
% isonormals(X/d,Y/d,Z/d,T,pgly1)
% pgly1.FaceColor = 'inTerp';
% pgly1.EdgeColor = 'none';
% daspect([1 1 1])
% axis(volumebounds(X/d,Y/d,Z/d,U,V,W))
% grid on
% colormap(flipud(gray))
% camlight
% lighting gouraud
% 
% clear T U V W X Y Z
% xlabel X/d_j
% zlabel Z/d_j
% grid on
% axis equal tight
% view(0,0)
% camup([-1 0 0])
% 
