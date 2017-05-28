function errest = error_estL(z_grid,z_edges,z,no_panels,rho,w,dzdt)
% z_grid: computational points for estimates
% z_edges: edges of our panels
% no_panels: nbr of panels
% rho: density from double layer
% w: quadrature weights
% dzdt: first derivative of boundary disc. points
panel_order = 16;

% Compute panel quantities for estimates
[L, numerator] = deal(zeros(1,no_panels));
for panel_idx=1:no_panels
    p_idx = (panel_idx-1)*panel_order+(1:panel_order);
    numerator(panel_idx) = max(abs(rho(p_idx).*imag(dzdt(p_idx))));
    L(panel_idx) = sum(abs(w(p_idx).*dzdt(p_idx)));
end
% Compute error estimate on grid
errest = zeros(size(z_grid));
tic
parfor idx=1:numel(z_grid)
    zi = z_grid(idx);
    % find closest panel edge
    [~, panel_idx] = min(abs(zi-z_edges));
    tmp = 0;
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
        tmp = tmp + 2*pi*numerator(panel_idx) ...
            *abs(z0 + 1i*sqrt(1-z0^2))^(-2*panel_order-1);
    end
    errest(idx) = tmp;
end
end