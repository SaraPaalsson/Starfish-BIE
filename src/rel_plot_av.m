%% PLOTTA MED REL.PROX

er = abs(u_known-uspec)/norm(u_known,Inf);

M = length(relp);

er_av = zeros(M,1);
for j=1:M
   indie = (j-1)*N+1:N*j;
   er_av(j) = norm(u_known(indie)-uspec(indie))/norm(u_known);
end

hold on
loglog(relp,er_av,'-.')
grid on
title('Error')
xlabel('Relative proximity')
ylabel('Relative Euclidean error')



max(er)
find(er == max(er))




