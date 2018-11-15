figure(1); clf;
plot(zDrops,'k-'); 
hold on; 
axis equal

%%

figure(1);
plot(z(il),'*')

%%

figure(2); clf; 
plot(tz,'k.-','MarkerSize',20); hold on; axis equal;



figure(2)
plot(z(il),'*')



figure(3); clf;
plot(nzpan,'k.-','MarkerSize',20);
hold on
plot(nzpan,0*nzpan,'r-')


figure(3)
plot(nz,'*')