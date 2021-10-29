%This Matlab script can be used to generate the variances of the Fourier
%random coefficients in the Fourier plane-wave series expansion of a 2D
%channel in Eq.(43). The script is valid for isotropic channels only and is 
%based on the theoretical computation in Appendix IV.C (part I) of the article:
%
%A. Pizzo, T. L. Marzetta and L. Sanguinetti, "Spatially-Stationary Model
%for Holographic MIMO Small-Scale Fading," in IEEE Journal on Selected Areas
%in Communications, vol. 38, no. 9, pp. 1964-1979, Sept. 2020,
%doi: 10.1109/JSAC.2020.3000877.
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


clear;
close all;
clc;

%% Parameters
%array size in number of wavelenghts (must be integer)
Lx = 16;

%% Variances of Fourier random coefficients
%discrete wavenumber frequencies
l_vec = [-Lx:1:Lx-1]';

%compute Fourier variances (2*Lx vector)
variances = asin((l_vec+1)/Lx) - asin(l_vec/Lx);

%normalize variances to their maximum
variances_norm = variances/max(variances);

%plot the variances in dB within the support segment
figure;
plot(l_vec,10*log10(variances_norm));
xlabel('$\ell$','Interpreter','Latex');
ylabel('$\sigma^2_{\ell}$ (dB)','Interpreter','Latex');
xlim([-Lx Lx])
grid on; box on;
set(gca,'FontSize',20);

