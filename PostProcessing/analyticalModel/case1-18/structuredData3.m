function field = structuredData3(datafile, gridfile, X, Y,Z, fieldStr)
%input datafile (.gfs) and gridfile (.dat) along with X Y and Z from mesh grid
%output field 

%% more about the gridfile: can be saved with this 

%gridfile='cartgrid2.dat';
%disp('saving the 2d grid');
%[X,Y,Z] = meshgrid(x,y,z);
%loc=[X(:),Y(:),Z(:)];
%save(gridfile,'loc','-ASCII','-SINGLE');

disp('interpolating and loading')    
tic, 
ll=evalc(['!gfs2oogl3D -p ' gridfile ' -c ' fieldStr ' < ' datafile ]); 
lolo=textscan(ll,'%f %f %f %f \n'); 
toc
try
    field=reshape(lolo{4},size(X)); 
catch
    disp('reshape failed while getting structured data.........')
    x=lolo{1}; y=lolo{2}; z = lolo{3};
    clear ll lolo
    plot(X(:),Y(:),Z(:),'ko',x,y,z,'r.')
    legend('Expected Points', 'Actual Points')
end
