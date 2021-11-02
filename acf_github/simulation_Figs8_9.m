%This Matlab script can be used to generate samples of the autocorrelation
%function of a 3D isotropic channel over a planar or volumetric array. The 
%autocorrelation function obtained with the Fourier plane-wave series expansion 
%is compared to the closed-form Clarke's isotropic autocorrelation function.
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
%array size in number of wavelenghts (must be integer, choose a power of 2 for IFFT)
Lx = 16;
Ly = 16;

%z coordinate normalized to the wavelength (must be an integer)
%for a planar array set z=0 (the forst entry must be 0)
z_vec = [0,1/2];

%number of Monte Carlo simulations
numOfMC = 1e4;

%%% Variances of Fourier random coefficients and migration coefficients
%discrete wavenumber frequencies
l_vec = [-Lx:1:Lx-1]';
m_vec = [-Ly:1:Ly-1];

%compute Fourier variances and migration coefficients (2*Ly x 2*Lx matrix)
[variances,kappaz] = function_computeVar(Lx,Ly);

%normalize variances so that they sum up to 1 (small-scale fading has unit power)
var_lm = variances/sum(variances(:));

%plot the variances in dB within the support ellipse
figure;
[X,Y] = meshgrid(l_vec,m_vec);
surf(X,Y,10*log10(var_lm/max(var_lm(:))));
colormap(parula(40))
colorbar
xlabel('$\ell$','Interpreter','Latex');
ylabel('$m$','Interpreter','Latex');
xlim([l_vec(1) l_vec(end)])
ylim([m_vec(1) m_vec(end)])
zlabel('$\sigma^2_{\ell,m}$ (dB)','Interpreter','Latex');
grid on; box on;
view(0,90);
set(gca,'FontSize',20);

%%% Generate 3D isotropic small-scale fading field over a plane - Clarke's method
%number of spatial samples (IFFT length) (Nyquist sampling:N>=2*L)
Nx = Lx*16; %high enough resolution
Ny = Ly*16; %high enough resolution
N = Nx*Ny;

acf_approx = zeros(N*size(z_vec,2),1);
for iter=1:numOfMC
    
    %update status simulation
    disp([num2str(iter),' /', num2str(numOfMC)]);
    
    %%% Generate 3D isotropic small-scale fading field  over a plane - Fourier plane-wave method
    %generate two 2D white noise random lattice fields with unit variance
    wlm_plus = sqrt(0.5)*(randn(2*Ly,2*Lx)+1i*randn(2*Ly,2*Lx));
    wlm_minus = sqrt(0.5)*(randn(2*Ly,2*Lx)+1i*randn(2*Ly,2*Lx));
    
    %generate two 2D independent random lattice fields with computed
    %variances - 2*Lx x 2*Ly matrix
    Hlm_plus = sqrt(var_lm).*wlm_plus/sqrt(2);
    Hlm_minus = sqrt(var_lm).*wlm_minus/sqrt(2);
    
    hnj_approx = [];
    %Go through all z-planes
    for indz=1:size(z_vec,2)
        
        z = z_vec(indz);
        
        %generate migration filters at z-plane
        migrationlm_plus = exp(1i*2*pi*kappaz*z);
        migrationlm_minus = exp(-1i*2*pi*kappaz*z);
        
        %spatial convolution over z
        Hlm_z = Hlm_plus.*migrationlm_plus + Hlm_minus.*migrationlm_minus;
        
        %apply zero-padding - Nx x Ny matrix
        Hlm_z_zeropad_l = [zeros(2*Ly,Nx/2-Lx), Hlm_z, zeros(2*Ly,Nx/2-Lx)];
        Hlm_z_zeropad = [zeros(Ny/2-Ly,Nx); Hlm_z_zeropad_l; zeros(Ny/2-Ly,Nx)];
        
        %prepare Fourier coefficients for IFFT - Nx x Ny matrix
        %     H_z = [Hlm_z(Ly+1:2*Ly,Lx+1:2*Lx),  Hlm_z(Ly+1:2*Ly,1:Lx);...
        %         Hlm_z(1:Ly,Lx+1:2*Lx),  Hlm_z(1:Ly,1:Lx)];
        Hlm_z_IFFT_l = fftshift(Hlm_z_zeropad,2);
        Hlm_z_IFFT = fftshift(Hlm_z_IFFT_l,1);
        
        %generate isotropic channel samples over a plane - Nx x Ny matrix
        hnj_approx_z = ifft2(Hlm_z_IFFT)*N;
        
        %save channel vectors at different z columnwise - N x 1 vector per
        %z element
        hnj_approx = [hnj_approx; hnj_approx_z(:)];
        
    end
    
    %compute spatial autocorrelation function - firt row of autocorrelation matrix
    acf_approx = acf_approx + real(hnj_approx(1)*conj(hnj_approx))/numOfMC;
    
end


%Go through all z-planes
for indz=1:size(z_vec,2)
    
    %select the autocorrelation function related to z
    acf_approx_z = acf_approx(1+(indz-1)*N:N*indz);
    
    %rearrange autocorrelation function
    acf_approx_z_mat = reshape(acf_approx_z,Ny,Nx);
    acf_approx_z_mat_l = ifftshift(acf_approx_z_mat,2);
    acf_approx_z = ifftshift(acf_approx_z_mat_l,1);
    
    %spatial grid
    x_axis = linspace(-Lx/2,Lx/2,Nx);
    y_axis = linspace(-Ly/2,Ly/2,Ny)';
    [y,x] = meshgrid(x_axis,y_axis);
    
    %plot spatial autocorrelation function
    figure;
    surf(y,x,acf_approx_z,'FaceAlpha',1);
    colormap(parula(100))
    cl = caxis;
    colorbar
    xlabel('$x/\lambda$','Interpreter','Latex')
    ylabel('$y/\lambda$','Interpreter','Latex')
    zlabel('$c(x,y)$','Interpreter','Latex')
    xlim([-Ly/4 Ly/4])
    ylim([-Lx/4 Lx/4])
    zl=zlim;
    grid on; box on;
    set(gca,'fontsize', 18)
    view(0,90)
    
    %Create a 3D uniform grid (normalized to the wavelength)
    x_axis_1 = (0:1:Nx/2-1)*Lx/Nx;
    y_axis_1 = (0:1:Ny/2-1)'*Ly/Ny;
    z = z_vec(indz);
    %
    Grid_x = repmat(x_axis_1,[Ny/2 1]);
    Grid_y = repmat(y_axis_1,[1 Nx/2]);
    arrayGeometry = [Grid_x(:), Grid_y(:), z*ones(N/4,1)];
    
    %Clarke's spatial autocorrelation function - firt row of autocorrelation matrix
    d_vec = vecnorm(arrayGeometry,2,2);
    acf_closedform_1 = reshape(sinc(2*d_vec),Ny/2,Nx/2);
    acf_closedform_col = [flip(acf_closedform_1,1);acf_closedform_1];
    acf_closedform_z = [flip(acf_closedform_col,2),acf_closedform_col];
    
    figure;
    surf(y,x,acf_closedform_z,'FaceAlpha',1);
    colormap(parula(100))
    caxis(cl)
    colorbar
    xlabel('$x/\lambda$','Interpreter','Latex')
    ylabel('$y/\lambda$','Interpreter','Latex')
    zlabel('$c(x,y)$','Interpreter','Latex')
    xlim([-Ly/4 Ly/4])
    ylim([-Lx/4 Lx/4])
    zlim(zl)
    grid on; box on;
    set(gca,'fontsize', 18)
    view(0,90)
    %legend('Closed-form','Approximated')
    
end


