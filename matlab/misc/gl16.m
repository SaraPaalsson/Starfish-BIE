function [zDrops, wDrops, tpar,panels] = gl16(Npanels,tau)
% Divide domain parametrized by tau into Npanels panels and distribute 16
% GL nodes and weights on each
%   zDrops, GL nodes in complex plane
%   spar, GL nodes in parametrization [0,2pi]
%   wDrops, GL weights

panels = linspace(0,2*pi,Npanels+1)'; %divide parametrization into panels

zDrops = zeros(Npanels*16,1); %complex interface coord.
tpar = zeros(Npanels*16,1); %parametrization vector
wDrops = zeros(Npanels*16,1); %GL. weights

for i=1:Npanels
    [nodes,weights] = gaussleg(16,[panels(i) panels(i+1)]);
    tpar((i-1)*16+1:16*i) = nodes;
    zDrops((i-1)*16+1:16*i) = tau(nodes);
    wDrops((i-1)*16+1:16*i) = weights;
end

end