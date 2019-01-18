% Estimate error Stokes

close all;
clear;
clc

% Either increase n or decrease imag(z0)
% casetest = 'incn';
% casetest = 'decz0';
casetest = 'linez0';

% Panel parametrisation (here flat!)
gamma = @(t) t;
gammap = @(t) t.^0;

mu = @(t) 2*sin(2*pi*t);
% mu = @(t) ones(size(t));

% qth derivative of k_n(z):
knq = @(n, q, z) (-(2*n+1)./sqrt(z.^2-1)).^q .* 2*pi./(z+sqrt(z.^2-1)).^(2*n+1);

switch casetest
    case 'decz0'
        n = 16;
        % Evaluation point in domai
        al = linspace(-2,2,20)';
        bl = flipud(linspace(0.001,1.5,100)');
        [a,b] = meshgrid(al,bl);
        z0mat = a + 1i*b;
        b = z0mat(:);
    case 'linez0'
        n = 16;
        a = 0.01;
        bl = flipud(linspace(0.001,1.5,20)');
        b = a + 1i*bl;
    case 'incn'
        z0 = 0.1 + 0.05i;
        nmax  = ceil(25/abs(imag(z0))); % Roughly gets us to full precision
        b = round(linspace(2, nmax, 100));
        b = b(1:2:end);
end


I1 = zeros(size(b)); I1corr = I1; I1err = I1; I1est = I1;
I2 = I1; I2corr = I1; I2err = I1; I2est = I1;
I3 = I1; I3corr = I1; I3err = I1; I3est = I1;
Itot_corr = I1; Itot_corr1 = I1; Itot_corr2 = I1; Itot_corr3 = I1;
Itot_err = I1; Itot_err1 = I1; Itot_err2 = I1; Itot_err3 = I1;
Itot_est = I1; Itot_est1 = I1; Itot_est2 = I1; Itot_est3 = I1;
Itot = I1; Itot_1 = I1; Itot_2 = I1; Itot_3 = I1;
for j = 1:length(b)
    
    switch casetest
        case 'decz0'
            %             z0 = a + 1i*b(j);
            z0 = b(j);
        case 'incn'
            n = b(j);
        case 'linez0'
            z0 = b(j);
    end
    z0_2 = abs(real(z0)) + 1i*abs(imag(z0));
    
    [tn,wn] = gaussleg(n,[-1 1]);
    zn = gamma(tn);
    znp = gammap(tn);
    mun = mu(zn);
    
    % First integrand
    f1 = @(t) (mu(t)./(gamma(t)-z0).*gammap(t));
    % Second integrand
    f2 = @(t) (mu(t)./((gamma(t)-z0).^2).*gammap(t));
    % Third integrand
    f3 = @(t) mu(t).*conj(gamma(t)-z0)./((gamma(t)-z0).^2).*gammap(t);
    
    % Compute first integral
    I1(j) = sum(mun.*wn.*(1./(zn-z0)).*znp);
    I1corr(j) = integral(f1, -1, 1, 'AbsTol', eps, 'RelTol', eps);
    I1err(j) = abs(I1corr(j)-I1(j));
    % Compute second integral
    I2(j) = sum(mun.*wn.*(1./((zn-z0).^2)).*znp);
    I2corr(j) = integral(f2, -1, 1, 'AbsTol', eps, 'RelTol', eps);
    I2err(j) = abs(I2corr(j)-I2(j));
    % Compute third integral
    I3(j) = sum(mun.*wn.*conj(zn-z0)./((zn-z0).^2).*znp);
    I3corr(j) = integral(f3,-1,1,'AbsTol',eps,'RelTol',eps);
    I3err(j) = abs(I3corr(j)-I3(j));
    
    %     mum = max(abs(mun));
    mum = abs(mu(z0));
    % First estimate
    I1est(j) = abs(mum.*knq(n,0,z0_2));
    % Second estimate
    I2est(j) = abs(mum.*knq(n,1,z0_2)/gammap(z0));
    % Third estimate
    I3est(j) = abs(mum.*knq(n,1,z0_2)*(conj(gamma(conj(z0)))-conj(gamma(z0)))/gammap(z0));
    
    % Total integrand
    nt = @(t) -1i*gammap(t)./abs(gammap(t));
    ftot1 = @(t) f1(t);
    ftot2 = @(t) f2(t).*nt(t);
    ftot3 = @(t) f3(t);
    
    ntn = nt(tn);
    
    Itot_corr1(j) = integral(ftot1,-1,1,'AbsTol',eps,'RelTol',eps);
    Itot_corr2(j) = integral(ftot2,-1,1,'AbsTol',eps,'RelTol',eps);
    Itot_corr3(j) = integral(ftot3,-1,1,'AbsTol',eps,'RelTol',eps);
    Itot_corr(j) =  Itot_corr1(j) + Itot_corr2(j) + Itot_corr3(j);
    
    Itot_1(j) = sum(wn.*mun.*znp./(zn-z0));
    Itot_2(j) = sum(wn.*mun.*ntn.*znp./((zn-z0).^2));
    Itot_3(j) = sum(wn.*mun.*znp.*conj(zn-z0)./((zn-z0).^2));
    Itot(j) =  Itot_1(j) + Itot_2(j) + Itot_3(j);
    Itot_err(j) = abs(Itot(j)-Itot_corr(j));
    Itot_err1(j) = abs(Itot_1(j) - Itot_corr1(j));
    Itot_err2(j) = abs(Itot_2(j) - Itot_corr2(j));
    Itot_err3(j) = abs(Itot_3(j) - Itot_corr3(j));
    
    mum1 = mu(z0).*gammap(z0);
%     mum1 = max(abs(mun.*gammap(zn)));
%     mum1 = max(abs(mun))*max(abs(gammap(zn)));
    Itot_est1(j) = (knq(n,0,z0_2).*mum1./gammap(z0));
    mum2 = mu(z0).*gammap(z0).*(1i*conj(gammap(conj(z0)))).^2./abs(gammap(z0));
%     mum2 = max(abs(mun)).*max(abs(gammap(zn))).*max(abs(1i.*conj(gammap(conj(zn))).^2))./max(abs(gammap(zn)));
    Itot_est2(j) = (mum2.*knq(n,1,z0_2)/gammap(z0));
    mum3 = mu(z0).*gammap(z0);
%     mum3 = max(abs(mun.*gammap(zn)));
    Itot_est3(j) = (mum3.*knq(n,1,z0_2)*(conj(gamma(conj(z0)))-conj(gamma(z0)))/(gammap(z0)^2));
    Itot_est(j) = abs(Itot_est1(j)) + abs(Itot_est2(j)) + abs(Itot_est3(j));
    Itot_est1(j) = abs(Itot_est1(j));
    Itot_est2(j) = abs(Itot_est2(j));
    Itot_est3(j) = abs(Itot_est3(j));
    
    %     Itot_corr(j) = 1/(pi*1i)*integral(ftot1,-1,1,'AbsTol',eps,'RelTol',eps) - ...
    %         1/(2*pi)*conj(integral(ftot2,-1,1,'AbsTol',eps,'RelTol',eps)) - ...
    %         1/(2*pi)*conj(integral(ftot3,-1,1,'AbsTol',eps,'RelTol',eps));
    %     Itot(j) = 1/(pi*1i)*sum(wn.*(znp./(zn-z0)).*mun) - ...
    %        1/(2*pi)*conj(sum(wn.*mun.*conj(ntn).^2.*znp./(zn-z0))) - ...
    %        1/(2*pi)*conj(sum(wn.*mun.*conj(zn-z0).*znp./((zn-z0).^2)));
    %     Itot_err(j) = abs(Itot(j)-Itot_corr(j));
    
    %     mum1 = (mu(z0).*gammap(z0));
    %     mum2 = (mu(z0)*conj(nt(z0)).^2.*gammap(z0));
    %     mum3 = (mu(z0).*gammap(z0));
    %     Itot1 = abs(1/(1i*pi)*imag(knq(n,0,z0_2))*mum1./gammap(z0));
    %     Itot2 = abs(1/(2*pi)*conj(knq(n,1,z0_2)/gammap(z0)*mum2./gammap(z0)));
    %     Itot3 = abs(1/(2*pi)*conj(knq(n,1,z0_2)*(conj(gamma(conj(z0)))-conj(gamma(z0)))/gammap(z0)^2*mum3));
    %     Itot_est(j) =  Itot1 + Itot2 + Itot3;
    
end

switch casetest
    case 'incn'
        figure(1);
        clf();
        h1 = semilogy(b, I1err, '.', 'MarkerSize',15,'DisplayName', ['$f_1 = \frac{1}{\tau-z_0}$']);
        hold on
        semilogy(b, I1est, 'k--', 'DisplayName', 'Est1');%,'Color',h1.Color)
        
        h2 = semilogy(b, I2err, '.', 'MarkerSize',15,'DisplayName', ['$f_2 = \frac{1}{(\tau-z_0)^2}$']);
        semilogy(b, I2est, 'k:', 'DisplayName', 'Est2');%,'Color',h2.Color)
        
        h3 = semilogy(b, I3err, '.', 'MarkerSize',15,'DisplayName', ['$f_3 = \frac{\bar{\tau}-\bar{z_0}}{(\tau-z_0)^2}$']);
        semilogy(b, I3est, 'k-.', 'DisplayName', 'Est3');%,'Color',h3.Color)
        ylabel('Error')
        hl = legend('toggle'); set(hl,'interpreter','latex')
        grid on
        % ylim([eps(), 10])
        xlabel('n')
        title(['$z_0 = $ '  num2str(z0)])
        
        figure(2); clf;
        h1 = semilogy(b, Itot_err, '.', 'MarkerSize',15,'DisplayName', ['$f_{tot} = \frac{1}{\tau-z_0}$']);
        hold on
        semilogy(b,Itot_est,'k:','DisplayName','Est tot')
        hl = legend('toggle'); set(hl,'interpreter','latex')
        xlabel('n')

        
    case 'linez0'
        I1err = Itot_err1; I1est =Itot_est1;
        I2err = Itot_err2; I2est = Itot_est2;
        I3err = Itot_err3; I3est = Itot_est3;
        
        figure(1);
        clf();
        subplot(121);
        h1 = semilogy(imag(b), I1err, '.', 'MarkerSize',15,'DisplayName', ['$f_1 = \frac{1}{\tau-z_0}$']);
        hold on
        semilogy(imag(b), I1est, 'k--', 'DisplayName', 'Est1');%,'Color',h1.Color)
        
        h2 = semilogy(imag(b), I2err, '.', 'MarkerSize',15,'DisplayName', ['$f_2 = \frac{1}{(\tau-z_0)^2}$']);
        semilogy(imag(b), I2est, 'k:', 'DisplayName', 'Est2');%,'Color',h2.Color)
        
        h3 = semilogy(imag(b), I3err, '.', 'MarkerSize',15,'DisplayName', ['$f_3 = \frac{\bar{\tau}-\bar{z_0}}{(\tau-z_0)^2}$']);
        semilogy(imag(b), I3est, 'k-.', 'DisplayName', 'Est3');%,'Color',h3.Color)
        
        ylabel('Error')
        hl = legend('toggle'); set(hl,'interpreter','latex')
        grid on
        ylim([eps(), 10])
        xlabel('n')
        title(['$z_0 = $ '  num2str(z0)])
        
        subplot(122)
        h4 = semilogy(imag(b), Itot_err, '.', 'MarkerSize',15,'DisplayName', ['$f_{tot}$']); hold on
        semilogy(imag(b),Itot_est,'k-','DisplayName','Est tot');
        ylabel('Error')
        hl = legend('toggle'); set(hl,'interpreter','latex')
        grid on
        ylim([eps(), 10])
        xlabel('n')
        title(['$z_0 = $ '  num2str(z0)])
        
        
    case 'decz0'
        levels = -15:1:0;
        
        figure(1); clf;
        %         I1err_mat = reshape(I1err,size(z0mat));
        I1err_mat = reshape(Itot_err1,size(z0mat));
        contourf(real(z0mat),imag(z0mat),log10(I1err_mat),levels,'EdgeColor','none')
        %         I1est_mat = reshape(I1est,size(z0mat));
        I1est_mat = reshape(Itot_est1,size(z0mat));
        hold on
        contour(real(z0mat),imag(z0mat),log10(I1est_mat),levels,'k','linewidth',3)
        title('Error + estimate for I1')
        xlabel('$t$, $\Re(z)$'); ylabel('$\Im(z)$')
        cbar = colorbar;
        caxis([-15 0])
        
        figure(2); clf;
        %         I2err_mat = reshape(I2err,size(z0mat));
        I2err_mat = reshape(Itot_err2,size(z0mat));
        contourf(real(z0mat),imag(z0mat),log10(I2err_mat),levels,'EdgeColor','none')
        %         I2est_mat = reshape(I2est,size(z0mat));
        I2est_mat = reshape(Itot_est2,size(z0mat));
        hold on
        contour(real(z0mat),imag(z0mat),log10(I2est_mat),levels,'k','linewidth',3)
        title('Error + estimate for I2')
        xlabel('$t$, $\Re(z)$'); ylabel('$\Im(z)$')
        cbar = colorbar;
        caxis([-15 0])
        
        figure(3); clf;
        %         I3err_mat = reshape(I3err,size(z0mat));
        I3err_mat = reshape(Itot_err3,size(z0mat));
        contourf(real(z0mat),imag(z0mat),log10(I3err_mat),levels,'EdgeColor','none')
        %         I3est_mat = reshape(I3est,size(z0mat));
        I3est_mat = reshape(Itot_est3,size(z0mat));
        hold on
        contour(real(z0mat),imag(z0mat),log10(I3est_mat),levels,'k','linewidth',3)
        title('Error + estimate for I3')
        xlabel('$t$, $\Re(z)$'); ylabel('$\Im(z)$')
        cbar = colorbar;
        caxis([-15 0])
        
        figure(4); clf;
        Itoterr_mat = reshape(Itot_err,size(z0mat));
        Itotest_mat = reshape(Itot_est,size(z0mat));
        contourf(real(z0mat),imag(z0mat),log10(Itoterr_mat),levels,'EdgeColor','none')
        hold on
        contour(real(z0mat),imag(z0mat),log10(Itotest_mat),levels,'k','linewidth',3)
        title('Error + estimate for Itot')
        xlabel('$t$, $\Re(z)$'); ylabel('$\Im(z)$')
        cbar = colorbar;
        caxis([-15 0])
        
        
end
