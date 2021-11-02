function [variances, kappaz] = function_computeVar(Lx,Ly)

%Function Inputs
%Lx:    array size over x-axis normalized to the wavelength (must be integer)
%Ly:    array size over y-axis normalized to the wavelength (must be integer)

%Function Outputs
%variances:      variances of the Fourier random coefficients
%kappaz:         exponential term of migration filter

%1st quadrant
l_vec = [0:1:Lx-1]';
m_vec = [0:1:Ly-1]';

%Initialization
var_1quadrant = zeros(Ly,Lx);
kappaz_1quadrant = zeros(Ly,Lx);
%Wavenumber support disk discretization
for mind=1:Ly
    m = m_vec(mind);
    for lind=1:Lx
        l = l_vec(lind);
        %parametrize the angle k_phi
        k_phi1 = atan(m*Lx/(l+1)/Ly);
        k_phi2 = min(atan(m*Lx/l/Ly),atan((m+1)*Lx/(l+1)/Ly));
        k_phi3 = max(atan(m*Lx/l/Ly),atan((m+1)*Lx/(l+1)/Ly));
        k_phi4 = atan((m+1)*Lx/l/Ly);
        %compute the polar integrals
        var_lm = function_computeInt(l,m,Lx,Ly,k_phi1,k_phi2,k_phi3,k_phi4);
        %returns the standard deviations for each (l,m) in the 1st quadrant
        var_1quadrant(Ly-mind+1,lind) = var_lm;
        %Migration filter coefficients
        kappaz_1quadrant(Ly-mind+1,lind) = real(sqrt(1-(l/Lx)^2-(m/Ly)^2));
    end
end

%exploit the rotational invariance of the integrand to extend the result to
%the other 3 quadrants
%2nd quadrant (symmetric wrt the ky axis)
var_2quadrant = fliplr(var_1quadrant);
kappaz_2quadrant = fliplr(kappaz_1quadrant);

%3rd quadrant (symmetric wrt the ky and kx axis)
var_3quadrant = flipud(fliplr(var_1quadrant));
kappaz_3quadrant = flipud(fliplr(kappaz_1quadrant));

%4th quadrant (symmetric wrt the kx axis)
var_4quadrant = flipud(var_1quadrant);
kappaz_4quadrant = flipud(kappaz_1quadrant);

%standard deviations for each (l,p)
variances = [var_2quadrant, var_1quadrant; var_3quadrant, var_4quadrant];
kappaz = [kappaz_2quadrant, kappaz_1quadrant; kappaz_3quadrant, kappaz_4quadrant];

%Discard possible negative values due to numerical computation
variances(variances<0) = 0;

end

