function dc = nakaDerivative( c, rmax, n, c50 )
c = c+0.25;
dc = ((n.*rmax.*c.^(n-1))./(c.^n+c.^c50)) - (rmax.*c.^n.*(n.*c.^(n-1)+c.^(c50-1).*c50))./((c.^n+c.^c50).^2);
