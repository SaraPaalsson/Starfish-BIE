function errest = error_estL_stokes(z_grid,z_edges,z,no_panels,rho1,rho2,rho3,w,dzdt)
% z_grid: computational points for estimates
% z_edges: edges of our panels
% no_panels: nbr of panels
% rho1: f(z) for first integral
% rho2: f(z) for second integral
% rho3: f(z) for third integral
% w: quadrature weights
% dzdt: first derivative of boundary disc. points
panel_order = 16;

% Compute panel quantities for estimates
[L, numerator1] = deal(zeros(1,no_panels));
% numerator2 = numerator1; numerator3 = numerator1;
for panel_idx=1:no_panels
    p_idx = (panel_idx-1)*panel_order+(1:panel_order);
%     numerator1(panel_idx) = max(abs(rho1(p_idx).*imag(dzdt(p_idx))));
%     numerator2(panel_idx) = max(abs(rho2(p_idx).*(dzdt(p_idx))));
% %     numerator3(panel_idx) = max(abs(rho3(p_idx).*(dzdt(p_idx))));
    L(panel_idx) = sum(abs(w(p_idx).*dzdt(p_idx)));
end
% Compute error estimate on grid
errest = zeros(size(z_grid));
tic
parfor idx=1:numel(z_grid)
    zi = z_grid(idx);
    % find closest panel edge
    [~, panel_idx] = min(abs(zi-z_edges));
    tmp1 = 0; tmp2 = 0; tmp3 = 0;
    % Go through panels connected to edge and sum error
    for panel_idx=[-1 0]+panel_idx
        panel_idx = mod(panel_idx-1,no_panels)+1;
        za = z_edges(panel_idx);
        zb = z_edges(panel_idx+1);
        p_idx = (panel_idx-1)*panel_order+(1:panel_order);
        zp = z(p_idx);
        % find distance to closest point on panel
        d = min(abs(zi-zp));
        % estimate projection on panel
        n1 = real(zb-za); n2 = imag(zb-za);
        v1 = real(zi-za); v2 = imag(zi-za);
        z0_r = -1+2*(n1*v1+n2*v2)/(n1*n1+n2*n2);
        % compute error estimate (abs->imag shows high freq oscillations)
        z0 = z0_r + 1i*d*2/L(panel_idx);
        
        
        % qth derivative of k_n(z):
        knq = @(n, q, z) (-(2*n+1)./sqrt(z.^2-1)).^q .* 2*pi./(z+sqrt(z.^2-1)).^(2*n+1);

        
%         mumax3 = max(abs());
%         numerator3 = max(abs(rho3(p_idx).*(dzdt(p_idx)))*mumax3);

        num1 = max(abs(rho1(p_idx)));
        num2 = max(abs(rho2(p_idx)));
        num3 = max(abs(rho3(p_idx)));
        
        tmp1 = tmp1 + num1*knq(panel_order,0,z0);
        tmp2 = tmp2 + num2*knq(panel_order,0,z0);
%         tmp3 = tmp3 + num3*knq(panel_order,1,z0)*(z0-conj(z0));
        tmp3 = 0;

% %         tmp1 = tmp1 + 2*pi*numerator1(panel_idx) ...
% %             *abs(z0 + 1i*sqrt(1-z0^2))^(-2*panel_order-1);
% %         
% %         tmp2 = tmp2 + 2*pi*numerator2(panel_idx) ...
% %             *abs(z0 + 1i*sqrt(1-z0^2))^(-2*panel_order-1);
%         
% %         tmp3 = tmp3 + 2*pi*numerator3 ...
% %             *abs(z0 + 1i*sqrt(1-z0^2))^(-2*panel_order-1);

%         tmp3 = tmp3 + numerator3*abs(pi./((z0 + sqrt(z0^2-1))^(2*panel_order+1))...
%         *((2*panel_order+1)^2/(sqrt(z0^2-1)^2) + z0*(2*panel_order+1)/(sqrt(z0^2-1)^3)));
    end
    errest(idx) = abs(tmp1) + abs(tmp2) + abs(tmp3);
end
end