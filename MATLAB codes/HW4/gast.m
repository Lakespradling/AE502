%--------------------------------------------------------------------------
%
% gast: Greenwich Apparent Sidereal Time
%
% Inputs:
%   Mjd_UT1   Modified Julian Date UT1
%   Mjd_TT    Modified Julian Date TT
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2022/06/16   Meysam Mahooti
%
%--------------------------------------------------------------------------
function gstime = gast(Mjd_UT1,Mjd_TT)

global const

gstime = mod(gmst(Mjd_UT1) + EqnEquinox(Mjd_TT), const.pi2);

