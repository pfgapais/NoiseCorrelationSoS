% Simulation of basic Sum-of-squares reconstruction
% Computation of SNR based on Psy

close all;  % close all figures
clear;
clc;
load psy_mat % Noise covariance matrix (N x N matrix)
load csm % coil sensitivity maps: B1- field maps (nb_voxels x nb_coils matrix)
csm=csm.*1000;%Scale signal to be above 100

% User inputs----------------------------------------------------------------------------
nb_port = 4;  % number of channels
nb_sample = 1e3;   % number of samples (time domain)
vn = 5; % nominal noise RMS voltage in I and Q
dvn = 50; % noise voltage dispersion (%)
cmax = 0.5; % maximum coefficient for noise linear combination (0 to 1)
% End of user inputs---------------------------------------------------------------------

% Define noises--------------------------------------------------------------
vn_vec = vn * (1 + rand(1, nb_port) * 2 * dvn / 100 - dvn / 100);% Draw noise amplitudes
c_mat = rand(nb_port, nb_port) * 2 * cmax - cmax;% Draw noise linear combination coefficients
c_mat = c_mat - diag(diag(c_mat)) + eye(nb_port);% Replace diagonal entries by 1
in_mat = randn(nb_sample, nb_port);% I component of noise (gaussian)
qn_mat = randn(nb_sample, nb_port);% Q component of noise (gaussian)
n0_mat = in_mat + 1i * qn_mat;% Create complex-valued noise
nc_mat = n0_mat * c_mat;% Linear combination of noise

nvar_mat = noise2covar(nc_mat); % compute linear combination covariance matrix
fnorm_vec = sqrt(2) ./ sqrt(diag(nvar_mat)); % noise normalization factor
nn_mat = zeros(nb_sample, nb_port); % initialize normalized noise
for k = 1:nb_port % noise normalization (noise has unit variance)
    nn_mat = fnorm_vec(k) .* nc_mat;    
end
n_mat = zeros(nb_sample, nb_port); % initialize scaled noise
for k = 1:nb_port % noise scaling with vn_vec
    n_mat = vn_vec(k) .* nn_mat;    
end

% Generate correlated noise from uncorrelated gaussian noise and noise covariance matrix
ncor_mat = n0_mat * chol(psy_mat) / sqrt(2); % generate noise data according to psy_mat 
psy_mat_new = noise2covar(ncor_mat); % recompute noise covariance matrix for checking
dpsy_mat = psy_mat - psy_mat_new; % difference matrix

% Initialize
nb_voxels=size(csm,1);
sn_mat = zeros(nb_sample, nb_port);
sos_reg = zeros(nb_sample, nb_port);
snr_reg = zeros(1, nb_voxels);

% Construct sn_mat for each voxel

for i = 1:nb_voxels
    % Add noise to signals for each voxel
    sn_mat = csm(i, :) + ncor_mat;
    % Compute SOS for the constructed sn_mat
    sos_reg = sos(sn_mat);
    % Evaluate signal and noise
    signal = mean(sos_reg, 1); % mean over samples
    noise = std(sos_reg, 0, 1); % std deviation over samples
    % Compute SNR
    snr_reg(i) = signal / noise;
end

disp(psy_mat);
disp(psy_mat_new);
disp(psy_mat-psy_mat_new);
% Make sure that psy_mat and psy_mat_new are very close to each others
% (if not, increase nb_sample)

save psy_mat psy_mat;save psy_mat_new psy_mat_new;save snr_reg snr_reg;

% Definition of functions----------------------------------------------------------------
function vec_out = sos(mat_in)
    % sum-of-squares
    % input is a nb_sample x nb_port matrix
    % output is a column vector (nb_sample x 1)
    [nbs, nbp] = size(mat_in); % find nb_sample and nb_port automatically
    vec_out = zeros(nbs, 1); % initialize output
    vec_out = sum(mat_in.*conj(mat_in),2);
    vec_out = sqrt(real(vec_out) / nbp); % real() to remove residual imaginary part
end

function ccovar = noise2covar(noise_mat)
    % Compute noise covariance matrix from noise samples organized
    % as a M x nbport matrix (noise_mat), M being the number of samples
    [~, np] = size(noise_mat); % find number of ports
    ccovar = zeros(np, np); % initialize complex covariance matrix
    for i = 1:np
        for j = 1:np
            Xi = noise_mat(:, i);
            Xj = noise_mat(:, j);
            EXi = mean(Xi);
            EXj = mean(Xj);
            covar = mean((Xi - EXi) .* conj(Xj - EXj));
            ccovar(i, j) = covar; % fill covariance matrix
        end
    end
end
% End of definitions-----------------------------------------------------------------------

