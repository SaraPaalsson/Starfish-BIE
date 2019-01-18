function ludvigtest()
    % Curve parametrization
    %gamma = @(t) t + exp(1i*t);
    x = @(t) t + cos(t);
    y = @(t) sin(t);
    gamma = @(t) x(t) + 1i*y(t);
    gammap = @(t) 1 + 1i*exp(1i*t);    
    % Pole
    t0 = 0.1+0.05i;
    z0 = gamma(t0);
    % Test functions
    f1z = @(z) 1./(z-z0); % regular f_1
    f2z = @(z) 1./(z-z0).^2; % regular f_2
    f3z = @(z) conj(z-z0)./(z-z0).^2; % Sara's function
    % Integrand including parametrization
    f1 = @(t) f1z(gamma(t)) .* gammap(t);
    f2 = @(t) f2z(gamma(t)) .* gammap(t);
    f3 = @(t) f3z(gamma(t)) .* gammap(t);
    % qth derivative of k_n(z):
    knq = @(n, q, z) (-(2*n+1)./sqrt(z.^2-1)).^q .* 2*pi./(z+sqrt(z.^2-1)).^(2*n+1);   
    % Compute quadrature errors
    nmax = ceil(25/abs(imag(t0))); % Roughly gets us to full precision
    nlist = round(linspace(2, nmax, 50));
    R1 = quaderr(f1, nlist);
    R2 = quaderr(f2, nlist);
    R3 = quaderr(f3, nlist);
    % Compute estimates (see AQBX paper also)
    Est1 = abs( knq(nlist, 0, t0));
    Est2 = abs( knq(nlist, 1, t0) / gammap(t0));
    % Limit of numerator at t -> t0
    % factor = x(t0) - x(t0') - 1i*(y(t0) - y(t0')); % As derived
    factor = gamma(t0')' - gamma(t0)'; % Shorthand for it
    Est3 = Est2*abs(factor);
    % Plot it    
    clf();
    h1 = semilogy(nlist, R1, '.', 'DisplayName', ['f1=' func2str(f1z)]);
    hold on
    semilogy(nlist, Est1, '--', 'DisplayName', 'Est1', 'Color', h1.Color)
    h2 = semilogy(nlist, R2, '.', 'DisplayName', ['f2=' func2str(f2z)]);
    semilogy(nlist, Est2, '--', 'DisplayName', 'Est2', 'Color', h2.Color)    
    h3 = semilogy(nlist, R3, '.', 'DisplayName', ['f3=' func2str(f3z)]);    
    semilogy(nlist, Est3, '--', 'DisplayName', 'Est3', 'Color', h3.Color)    
    ylabel('Error')
    xlabel('n')
    legend('toggle')
    grid on
    ylim([eps(), 10])
end

function R = quaderr(f, nlist)
% Quadrature error computer
    I = integral(f, -1, 1, 'AbsTol', eps, 'RelTol', eps);
    Q = zeros(size(nlist));
    for i=1:numel(nlist)
        n = nlist(i);
        [t, w] = lgwt(n, -1, 1);
        Q(i) = sum(w.*f(t));
    end
    R = abs(Q-I);    
end


% ============= Pasted from lgwt.m: 
function [x,w]=lgwt(N,a,b)
% lgwt.m
%
% [x,w]=lgwt(N,a,b)
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
%
% This code is covered by the BSD license (see bottom of file).

    N=N-1;
    N1=N+1; N2=N+2;

    xu=linspace(-1,1,N1)';

    % Initial guess
    y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

    % Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);

    % Derivative of LGVM
    Lp=zeros(N1,N2);

    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method

    y0=2;

    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0))>eps
        
        
        L(:,1)=1;
        Lp(:,1)=0;
        
        L(:,2)=y;
        Lp(:,2)=1;
        
        for k=2:N1
            L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end
        
        Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
        
        y0=y;
        y=y0-L(:,N2)./Lp;
        
    end

    % Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;      

    % Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

    % Copyright (c) 2009, Greg von Winckel
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    %
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
end