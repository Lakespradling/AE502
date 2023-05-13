%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to 
%    - the Earth's harmonic gravity field, 
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   t           Time since epoch in [s]
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2022/10/25   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function [dY] = Accel(t, Y)

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

MJD_TDB = Mjday_TDB(Mjd_TT);
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE440(MJD_TDB);

% Acceleration due to harmonic gravity field
a = AccelHarmonic(Y(1:3), E, AuxParam.n, AuxParam.m);

% Luni-solar perturbations
if (AuxParam.sun)
    a = a + AccelPointMass(Y(1:3),r_Sun,const.GM_Sun);
end

if (AuxParam.moon)
    a = a + AccelPointMass(Y(1:3),r_Moon,const.GM_Moon);
end

% Planetary perturbations
if (AuxParam.planets)
    a = a + AccelPointMass(Y(1:3),r_Mercury,const.GM_Mercury);
    a = a + AccelPointMass(Y(1:3),r_Venus,const.GM_Venus);
    a = a + AccelPointMass(Y(1:3),r_Mars,const.GM_Mars);
    a = a + AccelPointMass(Y(1:3),r_Jupiter,const.GM_Jupiter);
    a = a + AccelPointMass(Y(1:3),r_Saturn,const.GM_Saturn);
    a = a + AccelPointMass(Y(1:3),r_Uranus,const.GM_Uranus);    
    a = a + AccelPointMass(Y(1:3),r_Neptune,const.GM_Neptune);
    a = a + AccelPointMass(Y(1:3),r_Pluto,const.GM_Pluto);
end

dY = [Y(4:6);a];

