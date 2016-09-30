er = abs(u_known-u)/norm(u_known,Inf);
Er = reshape(log10(er),size(Rplot));
er_spec = abs(u_known-uspec)/norm(u_known,Inf);
Erspec = reshape(log10(er_spec),size(Rplot));

error = {reshape(er,size(Rplot))  reshape(er_spec,size(Rplot))};

spec_str1 = sprintf('Quadrature error, %d panels, normal quadrature',Npanels);
spec_str2 = sprintf('Quadrature error, %d panels, special quadrature',Npanels);
titstr = {spec_str1 spec_str2};


for i=1:2
   figure(i)
   pcolor(real(zplot),imag(zplot),log10(error{i}))
   shading flat
   % shading interp
   colormap parula
   colorbar
   caxis([-18 -2])
   ylabel('10log error normal G.-L.')
   set(gca,'yaxislocation','right');
   set(gca,'YTicklabel',[])
   set(gca,'XTicklabel',[])
    title(titstr{i})
end



% figure(1)
% % subplot(212)
% pcolor(real(zplot),imag(zplot),Er)
% shading flat
% % shading interp
% colormap parula
% % colormap (1-gray)
% colorbar
% caxis([-18 -2])
% % title('10log error normal GL')
% ylabel('10log error normal G.-L.')
% set(gca,'yaxislocation','right');
% % axis off
% set(gca,'YTicklabel',[])
% set(gca,'XTicklabel',[])
% % 
% % %%
% % % close all
% % 
% 
% figure(2)
% pcolor(real(zplot),imag(zplot),Erspec)
% shading flat
% % shading interp
% colormap parula
% % colormap (1-gray)
% colorbar
% % title('10log error special quad')
% % caxis([-18 -2])
% ylabel('10log error special quad')
% set(gca,'yaxislocation','right');
% % axis off
% set(gca,'YTicklabel',[])
% set(gca,'XTicklabel',[])
