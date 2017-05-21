clc
clear
close
d = 0.005;
g = 9.8067;
Ldomain = 10*d;
error = zeros(218,1);
eta = zeros(218,1);

for caseNo = 1:1:218
    tic
    %% Flow parameters
    Fr = 2.5; % Froude number of the jet
    alpha = 30;  % angle in degree
    Bo = 6; % Bond number
    Re = 85; % Reynolds number of the flow
    uj =  Fr*sqrt(9.8067*d); % velocity in m/s
    %% domain parameters
    xc = -0.077*Ldomain; %% origin of the flow -- point where liquid jets collide
    xmin = xc; %% keep the minimum doamin location there
    xmax = 1*Ldomain; %% length of the first link

    nx = 2500;

    %% Fiting parameters
    o = 2; %offset
    %% rejection parameters
    xrej = 0;

    %% calculation part
    etavec = 1.85:0.001:10; %% space set for eta. If the error is not below, 10%; try a bigger set of eta prodictions.
    errorEta = zeros(length(etavec),1);
    for ieta = 1:1:length(etavec)
        etatemp = etavec(ieta);
        x = linspace(0, (xmax - xc)/d +o, nx);
        e = 1; % coefficient of restitution

        psi0 = atand(e*tand(alpha));
        b = etatemp*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
        a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
        zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
        z = zeros(length(x),1);

        for jint = 2:1:length(x)
            z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
        end
        xp = x - o;
        zp = z;
    %% Rejection Part
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

    %% corelation fit: loads the coefficients from regression analysis
        load('finalRegressAnalytical')
        xs = xp;
        Cp1 = exp(b1(1));
        n1p1 = b1(2);
        n2p1 = b1(3);
        n3p1 = b1(4);
        n4p1 = b1(5);
        p1 = Cp1*sind(alpha)^n1p1*Fr^n2p1*Bo^n3p1*(Re/Fr)^n4p1;

        Cp2 = exp(b2(1));
        n1p2 = b2(2);
        n2p2 = b2(3);
        n3p2 = b2(4);
        n4p2 = b2(5);
        p2 = Cp2*sind(alpha)^n1p2*Fr^n2p2*Bo^n3p2*(Re/Fr)^n4p2;

        Cp3 = exp(b3(1));
        n1p3 = b3(2);
        n2p3 = b3(3);
        n3p3 = b3(4);
        n4p3 = b3(5);
        p3 = Cp3*sind(alpha)^n1p3*Fr^n2p3*Bo^n3p3*(Re/Fr)^n4p3;

        Cp4 = exp(b4(1));
        n1p4 = b4(2);
        n2p4 = b4(3);
        n3p4 = b4(4);
        n4p4 = b4(5);
        p4 = Cp4*sind(alpha)^n1p4*Fr^n2p4*Bo^n3p4*(Re/Fr)^n4p4;

        zf = p1*xs.^3 - p2*xs.^2 + p3*xs + p4;

        errorTemp = abs((zf-zp)./zf);
        errorEta(ieta) = sum(errorTemp)/length(zf); % L1 error measure
        msg = sprintf('Completed %d of %d in caseNo %d',ieta,length(etavec),caseNo);
        disp(msg)
    end
    [error(caseNo),id] = min(errorEta);
    etaf = etavec(id);
    eta(caseNo) = etaf;
    %% check and plot: Once eta is found! Plot the results to check
    x = linspace(0, (xmax - xc)/d +o,nx);
    e = 1; % coefficient of restitution
    psi0 = atand(e*tand(alpha));
    b = etaf*sind(psi0); %% b = (A + sin(\alpha)) make A some fraction of sin(\alpha)... can be more than 2 too
    a = (uj^2*(cosd(alpha)^2 + e^2*sind(alpha)^2))/(2*g);
    zd = tan(asin(sind(psi0) + b*(1./sqrt(x*d/(a) + 1) -1 )));
    z = zeros(length(x),1);

    for jint = 2:1:length(x)
        z(jint) = trapz(x(1:jint)*d,zd(1:jint))/d;
    end
    xp = x - o;
    zp = z;
    %% Rejection Part
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

    %% corelation fit
    load('finalRegressAnalytical')
    xs = xp;
    Cp1 = exp(b1(1));
    n1p1 = b1(2);
    n2p1 = b1(3);
    n3p1 = b1(4);
    n4p1 = b1(5);
    p1 = Cp1*sind(alpha)^n1p1*Fr^n2p1*Bo^n3p1*(Re/Fr)^n4p1;

    Cp2 = exp(b2(1));
    n1p2 = b2(2);
    n2p2 = b2(3);
    n3p2 = b2(4);
    n4p2 = b2(5);
    p2 = Cp2*sind(alpha)^n1p2*Fr^n2p2*Bo^n3p2*(Re/Fr)^n4p2;

    Cp3 = exp(b3(1));
    n1p3 = b3(2);
    n2p3 = b3(3);
    n3p3 = b3(4);
    n4p3 = b3(5);
    p3 = Cp3*sind(alpha)^n1p3*Fr^n2p3*Bo^n3p3*(Re/Fr)^n4p3;

    Cp4 = exp(b4(1));
    n1p4 = b4(2);
    n2p4 = b4(3);
    n3p4 = b4(4);
    n4p4 = b4(5);
    p4 = Cp4*sind(alpha)^n1p4*Fr^n2p4*Bo^n3p4*(Re/Fr)^n4p4;

    zf = p1*xs.^3 - p2*xs.^2 + p3*xs + p4;
    %% Final plotting for checking
    % Plot fit with data.
    figure( 'Visible', 'off' );
    if ~exist('fitplots', 'dir')
    mkdir('fitplots');
    end
    plot(xs,zf,'k-','LineWidth',3)
    hold on
    plot(xs,-zf,'k-','LineWidth',3)
    plot(xp,zp,'g-',xp,-zp,'r-','markers',15)
    hold on
    camup([-1 0 0])
    %axis([0 15 -1.5 1.5])
    axis equal
    xlabel('X/d_j')
    ylabel('Y/d_j')
    hold on
    name = sprintf('fitplots/%d.png',caseNo);
    saveas(gcf,name);
    close
    msg = sprintf('Case %d of 218 done',caseNo);
    disp(msg)
    toc

end
% save eta variable: this will be used in regression
save('eta.mat','eta','error')
