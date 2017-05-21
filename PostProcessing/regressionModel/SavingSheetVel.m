%% finding out the sheet average velocity

%% initial formalities
clc
clear
close

%% global parameters
place = 'intermediate/sim0.250.gfs';
gridfile='cartgrid.dat';
savingFile = 'sheeetVelocity.mat';
Fr = 2.5;

%% Geometrical Parameters
d = 0.005;
uj = Fr*sqrt(9.80665*d);
qj = 2*uj*pi/4*d^2;

%% Domain parameters
alpha = 30;
Ldomain = 10*d;
xmin = -0.1*Ldomain;
%xc = 0;
xmax = 1.36*Ldomain;
xc = -0.07745*Ldomain + 0.5*d/sind(alpha);
zmin = -Ldomain/6;
zmax = Ldomain/6;
ymin = 0;
ymax = 0.0325*Ldomain;
nth = 181;
nr = 60;
ny = 1000;
y = linspace(ymin,ymax,ny);
qf = zeros(nth,1);
ur = zeros(nth,1);
uth = zeros(nth,1);
th = linspace(0,pi,nth);


%% travering in theta
for thc = 1:1:nth
    sum1 = 0; % u - velocity : q flux
    sum2 = 0; % v - veocity : q flux
    hr = 0; % thickness times r
    count = 0; % counter for non-zero tracer in r
    getrMax; % sets the value of the maximum radius to travel in the post processing doamin box
    r = linspace(0,rmax,nr);
    for rc = 1:1:nr
        x = xc + r(rc)*cos(th(thc));
        z = r(rc)*sin(th(thc));
        [Y,X,Z] = meshgrid(y,x,z);
        loc = [X(:) Y(:) Z(:)];
        save(gridfile,'loc','-ASCII','-SINGLE');
        clear loc
        tic
        display(sprintf('At x/(l) = %4.3f and z/(l) = %4.3f : angle of %6.3f',x/Ldomain,z/Ldomain, th(thc)/pi*180))
        T = structuredData(place,gridfile,Y,X,'T');
        %break;
        if nnz(T) == 0
            display(sprintf('Got no non zero tracer at x/(l) = %4.3f and z/(l) = %4.3f',x/Ldomain,z/Ldomain, th(thc)/pi*180))
            toc
            break;
        else
            display('All good');
            count = count + 1;
            U = structuredData(place,gridfile,Y,X,'U');
            W = structuredData(place,gridfile,Y,X,'W');
            sum1 = sum1 + sum(U.*T)/sum(T); % integral of udy from 0 to h divided by h
            sum2 = sum2 + sum(W.*T)/sum(T); % integral of vdy from 0 to h divided by h
            hr = hr + sum(y.*T)/sum(T)*r(rc); % thickness times r
            toc
        end % end if of tracer equal to zero through out
        %break;
    end % end r - direction travel

    if count == 0
        display('Mayday Mayday! Something terrible happened')
        ur = 0;
        uth = 0;
        qf = 1;
        break;
    end % this should never happen

    uavg = (sum1*cos(th(thc)) + sum2*sin(th(thc)))/count;
    hr = hr/count;
    ur(thc) = uavg;
    uth(thc) = (- sum1*sin(th(thc)) + sum2*cos(th(thc)))/count;
    qf(thc) = uavg.*hr/qj;
    %break;
end % end theta - direction travel
hrm = qf./ur;
th = linspace(0,pi,181)';
urf = createFit(th,ur,1);
uthf = createFit(th,uth,2);
u = sqrt(urf(th).^2 + uthf(th).^2);
uf = createFit1(th,u,3);
u = mean(uf(th))/uj;

save(savingFile,'u','th','uf','hrm','qf')
