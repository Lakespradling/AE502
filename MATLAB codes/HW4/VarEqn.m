%------------------------------------------------------------------------------
%
% VarEqn.m
%
% Purpose:
%   Computes the variational equations, i.e. the derivative of the state vector
%   and the state transition matrix
%
% Input:
%   t           Time since epoch in [s]
%   yPhi        (6+36)-dim vector comprising the state vector (y) and the
%               state transition matrix (Phi) in column wise storage order
%
% Output:
%   yPhip       Derivative of yPhi
% 
% Last modified:   2022/10/25   Meysam Mahooti
%
%------------------------------------------------------------------------------
function [yPhip] = VarEqn(t, yPhi)

global const AuxParam eopdata

Mjd_UTC = AuxParam.Mjd_UTC+t/86400;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
Mjd_TT = Mjd_UTC + TT_UTC/86400;

% % Transformation matrix
% P = PrecMatrix(const.MJD_J2000,Mjd_TT);
% N = NutMatrix(Mjd_TT);
% T = N * P;
% E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1,Mjd_TT) * T;
E = GHAMatrix(Mjd_UT1,Mjd_TT);

% State vector components
r = yPhi(1:3);
v = yPhi(4:6);
Phi = zeros(6);

% State transition matrix
for j=1:6
    Phi(:,j) = yPhi(6*j+1:6*j+6);
end

% Acceleration and gradient
a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
G = EarthHarmGrad(r, E, AuxParam.n, AuxParam.m);

% Time derivative of state transition matrix
yPhip = zeros(42,1);
dfdy = zeros(6);

for i=1:3
    for j=1:3
        dfdy(i,j) = 0.0;                 % dv/dr(i,j)
        dfdy(i+3,j) = G(i,j);            % da/dr(i,j)
        if ( i==j )
            dfdy(i,j+3) = 1;
        else
            dfdy(i,j+3) = 0;             % dv/dv(i,j)
        end
        dfdy(i+3,j+3) = 0.0;             % da/dv(i,j)
    end
end

Phip = dfdy*Phi;

% Derivative of combined state vector and state transition matrix
for i=1:3
    yPhip(i)   = v(i);                 % dr/dt(i)
    yPhip(i+3) = a(i);                 % dv/dt(i)
end

for i=1:6
    for j=1:6
        yPhip(6*j+i) = Phip(i,j);     % dPhi/dt(i,j)
    end
end
  
