clc
clear
close

d = 0.005;
Ldomain = 10*d;
place = 'intermediate/sim0.250.gfs'; %% intermediate gerris simulation file
alpha = 30; % angle of impingement
xmin = -0.1*Ldomain;
xc = -0.077*Ldomain + 0.5*d/sind(alpha);
xmax = 1.3375*Ldomain;

zmin = 0;
zmax = Ldomain/5;
nx = 1200;
nz = 180;
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
toc
X = (X-xc)/d;
Z = Z/d;
[c,h] = contour(X,Z,T,[0.5 0.5],'k-','LineWidth',2);
close
if c(2,1) < 3*nx/4
    xs = zeros(nx,1);
    zs = zeros(nx,1);
    for icol = 1:1:nx
        tempold = 0;
        for irow = 1:1:nz
            if T(irow,icol) > 0.5
                temp = Z(irow,icol);
                if temp > tempold
                    tempold = temp;
                end
            end
        end
        zs(icol) = tempold;
        xs(icol) = X(irow,icol);
    end
else
    xs = c(1,2:c(2,1)+1);
    zs = c(2,2:c(2,1)+1);
end
[hf,gof] = createFitHx(xs,zs);
save('regression.mat','xs','zs','hf')
if gof.rsquare < 0.965
    [c,h] = contour(X,Z,T,[0.5 0.5],'k-','LineWidth',2);
    hold on
    plot(xs,hf(xs),'r-','LineWidth',2)
    xlabel X/d_j
    ylabel Z/d_j
    grid on
    axis equal
    camup([-1 0 0])
    if ~exist('fitplots', 'dir')
    mkdir('fitplots');
    end
    name = sprintf('fitplots/fitError.png');
    saveas(gcf,name);
    close
else
    [c,h] = contour(X,Z,T,[0.5 0.5],'k-','LineWidth',2);
    hold on
    plot(xs,hf(xs),'r-','LineWidth',2)
    xlabel X/d_j
    ylabel Z/d_j
    grid on
    axis equal
    camup([-1 0 0])
    if ~exist('fitplots', 'dir')
    mkdir('fitplots');
    end
    name = sprintf('fitplots/fitOK.png');
    saveas(gcf,name);
    close
end
clear X Z T loc c h
