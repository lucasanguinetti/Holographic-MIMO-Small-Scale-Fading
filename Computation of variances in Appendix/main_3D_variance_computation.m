%This Matlab script can be used to generate the variances of the Fourier
%random coefficients in the Fourier plane-wave series expansion of a 3D
%channel in Eq.(39). The script is valid for isotropic channels only and is 
%based on the theoretical computation in Appendix IV.C (part II) of the article:
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
Ly = 16;

%% Variances of Fourier random coefficients
%discrete wavenumber frequencies
l_vec = [-Lx:1:Lx-1]';
m_vec = [-Ly:1:Ly-1];

%compute Fourier variances (2*Ly x 2*Lx matrix)
[variances,~] = function_computeVar(Lx,Ly);

%normalize variances to their maximum
variances_norm = variances/max(variances(:));

%plot the variances in dB within the support ellipse
figure;
[X,Y] = meshgrid(l_vec,m_vec);
surf(X,Y,10*log10(variances_norm));
colormap(parula(40))
colorbar
xlabel('$\ell$','Interpreter','Latex');
ylabel('$m$','Interpreter','Latex');
xlim([-Lx Lx])
ylim([-Ly Ly])
zlabel('$\sigma^2_{\ell,m}$ (dB)','Interpreter','Latex');
grid on; box on;
view(0,90);
set(gca,'FontSize',20);

