%This Matlab script can be used to generate samples of the autocorrelation
%function of a 2D isotropic channel over a linear array. The autocorrelation
%function obtained with the Fourier plane-wave series expansion is compared
%to the closed-form Jake's isotropic autocorrelation function.
%The script is based on the Fourier plane-wave series expansion introduced
%in in Section V of the article:
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

%%% Parameters
%array size in number of wavelenghts (must be an integer)
Lx = 16;

%y coordinate normalized to the wavelength (must be an integer)
%for a planar array set y=0 (the forst entry must be 0)
y_vec = [0,1];

%number of Monte Carlo simulations
numOfMC = 1e4;

%%% Variances of Fourier random coefficients and migration coefficients
%discrete wavenumber frequencies
l_vec = [-Lx:1:Lx-1];

%compute Fourier variances (2*Lx vector)
%normalize variances so that they sum up to 1 (small-scale fading has unit power)
var_l = (asin((l_vec+1)/Lx) - asin(l_vec/Lx))/pi;

%compute migration coefficients (2*Lx vector)
kappay = sqrt(1-(l_vec/Lx).^2);

%plot the variances in dB within the support segment
figure;
plot(l_vec,10*log10(var_l/max(var_l)),'LineWidth',2);
xlabel('$\ell$','Interpreter','Latex');
ylabel('$\sigma^2_{\ell}$ (dB)','Interpreter','Latex');
xlim([-Lx Lx])
grid on; box on;
set(gca,'FontSize',20);

%%% Generate 2D isotropic small-scale fading field over a line - Jake's method
%number of spatial samples (IFFT length) (Nyquist sampling:N>=2*L)
Nx = Lx*16; %high enough resolution

acf_approx = zeros(Nx*size(y_vec,2),1);
for iter=1:numOfMC
    
    %update status simulation
    disp([num2str(iter),' /', num2str(numOfMC)]);
    
    %%% Generate a 1D isotropic small-scale fading field over a line - Fourier plane-wave method
    %generate two independent 1D white noise random lattice fields with unit variance
    wl_plus = sqrt(0.5)*(randn(1,2*Lx)+1i*randn(1,2*Lx));
    wl_minus = sqrt(0.5)*(randn(1,2*Lx)+1i*randn(1,2*Lx));
    
    %generate two 1D independent random lattice fields with computed variances - 2*Lx vector
    Hl_plus = sqrt(var_l).*wl_plus/sqrt(2);
    Hl_minus = sqrt(var_l).*wl_minus/sqrt(2);
    
    hn_approx = [];
    %Go through all y-planes
    for indy=1:size(y_vec,2)
        
        y = y_vec(indy);
        
        %generate migration filters at y-plane
        migrationl_plus = exp(1i*2*pi*kappay*y);
        migrationl_minus = exp(-1i*2*pi*kappay*y);
        
        %spatial convolution over y
        Hl_y = Hl_plus.*migrationl_plus + Hl_minus.*migrationl_minus;
        
        %apply zero-padding - Nx vector
        Hl_zeropad = [zeros(1,Nx/2-Lx), Hl_y, zeros(1,Nx/2-Lx)];
        
        %prepare Fourier coefficients for IFFT
        %Hl_IFFT = [Hl(Lx+1:2*Lx),  Hl(1:Lx)];
        Hl_IFFT = fftshift(Hl_zeropad);
        
        %generate isotropic channel samples over a line
        hn_approx_y = ifft(Hl_IFFT)*Nx;
        
        %save channel vectors at different y - Nx vector per y element
        hn_approx = [hn_approx; hn_approx_y.'];
        
    end
    
    %compute spatial autocorrelation function - firt row of autocorrelation matrix
    acf_approx = acf_approx + real(hn_approx(1)*conj(hn_approx))/numOfMC;
    
end

%Go through all y-planes
for indy=1:size(y_vec,2)
    
    %select the autocorrelation function related to z
    acf_approx_y = acf_approx(1+(indy-1)*Nx:Nx*indy);
    
    %rearrange autocorrelation function
    acf_approx_y = reshape(acf_approx_y,1,Nx);
    
    %rearrange autocorrelation function
    acf_approx_y = ifftshift(acf_approx_y);
    
    %Create a 3D uniform grid (normalized to the wavelength)
    grid_x = (-Nx/2:1:Nx/2-1)*Lx/Nx;
    y = y_vec(indy);
    %
    arrayGeometry = [grid_x; y*ones(1,Nx)];
    
    %Clarke's spatial autocorrelation function - firt row of autocorrelation matrix
    acf_closedform_y = besselj(0,2*pi*vecnorm(arrayGeometry,2,1));
    
    %%% Plot spatial autocorrelation function
    figure;
    plot(grid_x,acf_closedform_y,'-k','LineWidth',2); hold on
    plot(grid_x,acf_approx_y,'-.r','LineWidth',2);
    xlabel('$x/\lambda$','Interpreter','Latex')
    ylabel('$c(x)$','Interpreter','Latex')
    xlim([-Lx/4 Lx/4])
    ylim([-.5 1.1])
    legend('Clarke model','Fourier model')
    grid on; box on;
    set(gca,'fontsize', 18)
        
end
