%% This takes in the data from a gerris intermediate file case1-18/sim0.250.gfs
%% the data is then converted on to a uniform structured grid from the non-uniform
%% octree grid. 
clc
clear
close
place = 'case1-18/sim0.250.gfs';
d = 0.005;

alpha = 30;
Ldomain = 10*d;
xmin = -0.1*Ldomain;
xc = -0.077*Ldomain + 0.5*d/sind(alpha);
xmax = 0.875*Ldomain;

zmin = -Ldomain/8;
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
T = structuredData(place, gridfile, X, Z, 'T'); % volume fraction (VOF) data converted for the uniform grid XYZ
U = structuredData(place, gridfile, X, Z, 'U');
W = structuredData(place, gridfile, X, Z, 'W');
UW = U;
WW = W;
for i = 1:1:size(T,1)
    for j = 1:1:size(T,2)
        if T(i,j) < 0.5
            UW(i,j) = 0;
            WW(i,j) = 0;
        end

     end
end
X = X-xc;
[c,h] = contour(X,Z,T,[0.5 0.5]);
close
%% finds position for interface
xs = c(1,2:c(2,1)+1)/d;
zs = c(2,2:c(2,1)+1)/d;
save('analytical.mat','xs','zs','X','Y','Z','UW','WW','d') %% Saves the necessary variables.
clear
close
