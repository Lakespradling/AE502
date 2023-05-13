%--------------------------------------------------------------------------
% 
%  AccelHarmonic: Computes the acceleration due to the harmonic gravity
%   			  field of the central body
% 
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [a] = AccelHarmonic(r, E, n_max, m_max)

global CS

gm    = 398600.4415e9; % [m^3/s^2]; GGM03C & GGM03S
r_ref = 6378.1363e3;   % Earth's radius [m]; GGM03C & GGM03S

% Body-fixed position 
r_bf = E * r;

% Auxiliary quantities
r_sqr =  dot(r_bf,r_bf);      % Square of distance
rho   =  r_ref*r_ref/r_sqr;
  
x0 = r_ref*r_bf(1)/r_sqr;     % Normalized
y0 = r_ref*r_bf(2)/r_sqr;     % coordinates
z0 = r_ref*r_bf(3)/r_sqr;

%
% Evaluate harmonic functions 
%   V_nm = (r_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
% and 
%   W_nm = (r_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
% up to degree and order n_max+1
%

% Calculate zonal terms V(n,0); set W(n,0)=0.0
V(1,1) = r_ref/sqrt(r_sqr);
W(1,1) = 0;
      
V(2,1) = z0 * V(1,1);
W(2,1) = 0;

for n=2:n_max+1
  V(n+1,1) = ( (2*n-1) * z0 * V(n,1) - (n-1) * rho * V(n-1,1) )/n;
  W(n+1,1) = 0;
end

% Calculate tesseral and sectorial terms 

for m=1:m_max+1
    
  % Calculate V(m,m) .. V(n_max+1,m)
  V(m+1,m+1) = (2*m-1) * ( x0*V(m,m) - y0*W(m,m) );
  W(m+1,m+1) = (2*m-1) * ( x0*W(m,m) + y0*V(m,m) );

  if m<=n_max
    V(m+2,m+1) = (2*m+1) * z0 * V(m+1,m+1);
    W(m+2,m+1) = (2*m+1) * z0 * W(m+1,m+1);
  end

  for n=m+2:n_max+1
    V(n+1,m+1) = ( (2*n-1)*z0*V(n,m+1) - (n+m-1)*rho*V(n-1,m+1) ) / (n-m);
    W(n+1,m+1) = ( (2*n-1)*z0*W(n,m+1) - (n+m-1)*rho*W(n-1,m+1) ) / (n-m);
  end

end

%
% Calculate accelerations ax,ay,az
%
ax = 0;
ay = 0;
az = 0;

for m=1:m_max+1
  for n=m:n_max+1
    if (m==1) 
      C = CS(n,1);   % = C_n,0
      ax = ax - C * V(n+1,2);
      ay = ay - C * W(n+1,2);
      az = az - (n)*C * V(n+1,1);
    
    else  
      C = CS(n,m);   % = C_n,m 
      S = CS(m-1,n); % = S_n,m 
      Fac = 0.5 * (n-m+1) * (n-m+2);
      ax = ax  + 0.5 * ( - C * V(n+1,m+1) - S * W(n+1,m+1) )...
              + Fac * ( + C * V(n+1,m-1) + S * W(n+1,m-1) );
      ay = ay  + 0.5 * ( - C * W(n+1,m+1) + S * V(n+1,m+1) )...
              + Fac * ( - C * W(n+1,m-1) + S * V(n+1,m-1) );
      az = az + (n-m+1) * ( - C * V(n+1,m)   - S * W(n+1,m) );
    end
  end
end

% Body-fixed acceleration
a_bf = (gm/(r_ref*r_ref))*[ax,ay,az]';

% Inertial acceleration 
a = E'*a_bf;

