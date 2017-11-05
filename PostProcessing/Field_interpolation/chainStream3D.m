  %% This file reads the intermediate file sim0.250.gfs and converts the structured non-uniform
%% data to uniform data in chain.mat
clc
clear
close all
place = 'intermediate/sim0.250.gfs';
folder = 'postprocess'; % output folder
opFolder = fullfile(cd, folder);
if ~exist(opFolder, 'dir')
mkdir(opFolder);
end
if ~exist(place,'file')
    disp('%file does not exist!..............................')
    disp(place)
else
    d = 0.005;
    Ldomain = 10*d;
    alpha = 30;
    xmin = -0.1*Ldomain;
    xc = -0.077*Ldomain + 0.5*d/sind(alpha);
    xmax = 1.1*Ldomain;
    zmin = -Ldomain/5;
    zmax = Ldomain/5;
    ymin = -0.15*Ldomain;
    ymax = 0.15*Ldomain;
    nx = 1200;
    ny = 150;
    nz = 200;
    x = linspace(xmin,xmax,nx);
    y = linspace(ymin,ymax,ny);
    z = linspace(zmin,zmax,nz);
    tic
    gridfile='cartgrid3.dat';
    disp('saving the 3d grid');
    [X,Y,Z] = meshgrid(x,y,z);
    loc=[X(:),Y(:),Z(:)];
    save(gridfile,'loc','-ASCII','-SINGLE');
    toc
    T = structuredData3(place, gridfile, X, Y,Z, 'T');
    U = structuredData3(place, gridfile, X, Y,Z, 'U');
    V = structuredData3(place, gridfile, X, Y,Z, 'V');
    W = structuredData3(place, gridfile, X, Y,Z, 'W');
    save('chain.mat','U','V','W','T','X','Y','Z')
end

fprintf('\n You job is finished.\n');
