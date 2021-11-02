function [y] = sinInt(a,k_phi)

%Function Inputs
%a:             auxiliary parameter specified below Eq. (72)
%k_phi:         angular wavenumber coordinate specified below Eq. (71)

%Function Outputs
%Y:             compute integral in Eq. (74)

y = 1/a*atan(cos(k_phi)/sqrt(a^2*sin(k_phi)^2-1)) - asin(cos(k_phi)/sqrt(1-1/a^2));

end



