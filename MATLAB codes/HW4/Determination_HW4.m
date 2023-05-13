%--------------------------------------------------------------------------
%  Initial Orbit Determination
%
%  Adapted from https://www.mathworks.com/matlabcentral/fileexchange/
%  55178-initial-orbit-determination-angles-only-method
%
% References:
%   O. Montenbruck, E. Gill, "Satellite Orbits: Models, Methods, and
%   Applications", Springer Verlag, Heidelberg, 2005.
%   
%   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
%   4th Edition, 2013.
%
%   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
%
% Last modified:   2023/05/12 Lake Spradling
%
%--------------------------------------------------------------------------
clc
clear all
format long g

global const CS AuxParam eopdata PC

SAT_Const
load DE440Coeff.mat
PC = DE440Coeff;
% Earth gravity field coefficients
GGM03C_unnormalized

% Model parameters
AuxParam = struct ('Mjd_UTC',0,'n',0,'m',0,'sun',0,'moon',0,'planets',0);

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
for i=1:23
    tline = fgetl(fid);
end
numrecsobs = str2num(tline(21:end));
tline = fgetl(fid);
for i=1:numrecsobs
    eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
end
for i=1:4
    tline = fgetl(fid);
end
numrecspred = str2num(tline(22:end));
tline = fgetl(fid);
for i=numrecsobs+1:numrecsobs+numrecspred
    eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
end
fclose(fid);

% Read observations
fid = fopen('OneWeb.txt','r');
nobs = 0; % number of observations
while ~feof(fid)
    nobs = nobs+1;
    tline = fgetl(fid);
    Y = str2num(tline(1:4));
    M = str2num(tline(6:7));
    D = str2num(tline(9:10));
    h = str2num(tline(13:14));
    m = str2num(tline(16:17));
    s = str2num(tline(19:24));
    az = str2num(tline(26:34));
    el = str2num(tline(37:end));
    obs(nobs,1) = Mjday(Y,M,D,h,m,s);
    % Subtract known biases
    obs(nobs,2) = const.Rad*(az-0.0081);
    obs(nobs,3) = const.Rad*(el-0.0045);
end
fclose(fid);

% Measurement uncertainty
sigma_az = 1e-5*const.Rad; % [rad]
sigma_el = 1e-5*const.Rad; % [rad]

% Champaign, IL station
lat = const.Rad*40.1164;    % [rad]
lon = const.Rad*(-88.2434);   % [rad]
alt = 233;                  % [m]

% Station's ECEF position vector
Rs = Position(lon, lat, alt)';
% find station's velocity vector
ome = [0; 0; const.omega_Earth];
Vs = cross(ome, Rs);

third = floor(nobs/2)+1;
second = floor(third/2)+1;

% % double-r-iteration method
% [rtod,vtod] = anglesdr(obs(1,2),obs(second,2),obs(third,2),obs(1,3),obs(second,3),obs(third,3),...
%                    obs(1,1),obs(second,1),obs(third,1),Rs,Rs,Rs);
% Y0_apr = [rtod;vtod];

% Gauss method
rtasc = zeros(nobs,1);
decl = zeros(nobs,1);
rsite = zeros(nobs,3);
for i=1:nobs
    [lst,gst] = lstime(lon,obs(i,1)+2400000.5);
    [rtasc(i),decl(i)] = azl2radc(obs(i,2),obs(i,3),lat,lst);
    RsVs = ECEF2ECI(obs(i,1),[Rs;Vs]');
    rsite(i,:) = 1e-3*RsVs(1:3);
end
[reci,veci] = anglesg(decl(1),decl(second),decl(third),rtasc(1),rtasc(second),...
                  rtasc(third),obs(1,1)+2400000.5,obs(second,1)+2400000.5,...
                  obs(third,1)+2400000.5,rsite(1,:),rsite(second,:),...
                  rsite(third,:),1e-3*const.R_Earth,1e-9*const.GM_Earth,86400);
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(second,1),'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT  = obs(second,1)+TT_UTC/86400;
T = (Mjd_TT-const.MJD_J2000)/36525;
[rtod,vtod,atod] = eci2tod(reci',veci',[0;0;0],T,dpsi,deps);
Y0_apr = 1e3*[rtod;vtod];

AuxParam.Mjd_UTC = obs(second,1);
AuxParam.n      = 20;
AuxParam.m      = 20;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;

options = rdpset('RelTol',1e-10,'AbsTol',1e-12);
[~,yout] = radau(@Accel,[0 -(obs(second,1)-obs(1,1))*86400],Y0_apr,options);
Y0 = yout(end,:)';

% Initialize the Least Squares' parameters
A = zeros((nobs-1)*2,6);
b = zeros((nobs-1)*2,1);
w = zeros((nobs-1)*2,(nobs-1)*2);
yPhi = zeros(42,1);
Phi  = zeros(6);

% Transformation from Greenwich meridian system to local tangent coordinates
LTC = LTCMatrix(lon,lat);

% Implementation of the Least Squares method
AuxParam.Mjd_UTC = obs(1,1);
for iterat = 1:10
    fprintf('\nIteration Nr. %d \n', iterat);
    fprintf('\nResiduals:\n');
    fprintf('    MjdUTC         Azim(deg)      Elev(deg)\n');
    
    for i=2:nobs
        % Time increment and propagation
        t = (obs(i,1)-obs(1,1))*86400; % Time since epoch [s]
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(i,1),'l');
        [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        Mjd_UT1 = obs(i,1) + UT1_UTC/86400;
        Mjd_TT = obs(i,1) + TT_UTC/86400;
        
        for ii=1:6
            yPhi(ii) = Y0(ii);
            for j=1:6  
                if (ii==j) 
                    yPhi(6*j+ii) = 1; 
                else
                    yPhi(6*j+ii) = 0;
                end
            end
        end
        
        [~,yout] = radau(@VarEqn,[0 t],yPhi,options);
        yPhi = yout(end,:)';
        
        % Extract state transition matrices        
        for j=1:6
            Phi(:,j) = yPhi(6*j+1:6*j+6);
        end

        [~,yout] = radau(@Accel,[0 t],Y0,options);
        Y = yout(end,:)';
        
        % Topocentric coordinates        
        U = GHAMatrix(Mjd_UT1,Mjd_TT);        % Earth rotation
        r = Y(1:3);
        s = LTC*(U*r-Rs);                     % Topocentric position [m]
        
        % Observations and partials
        [Azim, Elev, dAds, dEds] = AzElPa(s); % Azimuth, Elevation
        
        dAdY0 = [dAds*LTC*U,zeros(1,3)]*Phi;
        dEdY0 = [dEds*LTC*U,zeros(1,3)]*Phi;
        
        % Accumulate least-squares system
        A(2*(i-1)-1:2*(i-1),:) = [dAdY0;dEdY0];
        b(2*(i-1)-1:2*(i-1)) = [obs(i,2)-Azim ; obs(i,3)-Elev];
        w(2*(i-1)-1:2*(i-1),2*(i-1)-1:2*(i-1)) = diag([1/sigma_az^2,1/sigma_el^2]);
        
        fprintf(' %12f', obs(i,1));
        fprintf(' %14f', (obs(i,2)-Azim)*const.Deg);
        fprintf(' %14f\n', (obs(i,3)-Elev)*const.Deg);
    end
    
    % Solve least-squares system
    dY0 = inv(A'*w*A)*A'*w*b;
    
    fprintf('\n Correction: \n');
    fprintf(' Pos [m] \n');
    fprintf('%6.1f\n%6.1f\n%6.1f\n',dY0(1),dY0(2),dY0(3));
    fprintf('\n Vel [m/s] \n');
    fprintf('%6.1f\n%6.1f\n%6.1f\n',dY0(4),dY0(5),dY0(6));
    
    % Correct epoch state
    Y0 = Y0 + dY0;    
end

[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(1,1),'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT  = obs(1,1)+TT_UTC/86400;
T = (Mjd_TT-const.MJD_J2000)/36525;
[reci,veci,aeci] = tod2eci(Y0(1:3),Y0(4:6),[0;0;0],T,dpsi,deps);
Y0 = [reci;veci];
[p,a,ecc,incl,raan,argp,nu,m,arglat,truelon,lonper] = rv2coeh(1e-3*Y0(1:3),...
 1e-3*Y0(4:6),6378.1363,398600.4415);

