

er = abs(u_known-u)/norm(u_known,Inf);
Er = reshape(log10(er),size(Rplot));
figure()
% subplot(212)
pcolor(real(zplot),imag(zplot),Er)
shading interp
colormap parula
% colormap (1-gray)
colorbar
caxis([-18 -2])
% title('10log error normal GL')
ylabel('10log error normal G.-L.')
set(gca,'yaxislocation','right');
% axis off
set(gca,'YTicklabel',[])
set(gca,'XTicklabel',[])
% 
% %%
% % close all
% 
er_spec = abs(u_known-uspec)/norm(u_known,Inf);
Erspec = reshape(log10(er_spec),size(Rplot));
figure()
pcolor(real(zplot),imag(zplot),Erspec)
shading interp
colormap parula
% colormap (1-gray)
colorbar
% title('10log error special quad')
% caxis([-18 -2])
ylabel('10log error special quad')
set(gca,'yaxislocation','right');
% axis off
set(gca,'YTicklabel',[])
set(gca,'XTicklabel',[])
