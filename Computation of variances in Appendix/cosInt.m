function [y] = cosInt(a,k_phi)

%Function Inputs
%a:             auxiliary parameter specified below Eq. (72)
%k_phi:         angular wavenumber coordinate specified below Eq. (71)

%Function Outputs
%Y:             compute integral in Eq. (75)

y = -1/a*atan(sin(k_phi)/sqrt(a^2*cos(k_phi)^2-1)) + asin(sin(k_phi)/sqrt(1-1/a^2));

end
