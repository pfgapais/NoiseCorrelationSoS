% Simulation of basic Sum-of-squares reconstruction
% Computation of SNR based on Psy
% Validation of the procedure for 1 voxel

% Test with N channels: different signal and noise levels,
%                       different noise correlation levels.

close all;  % close all figures
clear;
clc;

% User inputs----------------------------------------------------------------------------
nb_port = 4;  % number of channels
nb_sample = 1e6;   % number of samples (time domain)
vs = 500; % nominal signal voltage
vn = 5; % nominal noise RMS voltage in I and Q
dvs = 20; % signal voltage dispersion (%)
dvn = 100; % noise voltage dispersion (%)
cmax = 0.5; % maximum coefficient for noise linear combination (0 to 1)
% End of user inputs---------------------------------------------------------------------


% Define signals and noises--------------------------------------------------------------
vs_vec = vs * (1 + rand(1, nb_port) * 2 * dvs / 100 - dvs / 100);% Draw signal amplitudes
ps_vec = rand(1, nb_port) * 2 * pi;% Draw signal phases
s_vec = vs_vec .* exp(1i * ps_vec);% Generate complex-valued signals
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

% Construct signal+noise data
sn_mat = zeros(nb_sample, nb_port); % initialize 
for k = 1:nb_port 
    sn_mat(:, k) = s_vec(k) + n_mat(:, k); % add noise to signals
end

% Compute direct sum-of-squares
sos_dir = sos(sn_mat); % compute SoS
signal = mean(sos_dir); % evaluate signal
noise = std(sos_dir); % evaluate noise
snr_dir = signal / noise; % compute snr
fprintf('SNR from direct computation (original noise) = %f\n', snr_dir);

% Plot sos_dir histogram
figure; % create graphic
histogram(sos_dir, 100, 'Normalization', 'pdf'); % plot histogram

% Generate correlated noise from uncorrelated gaussian noise and noise covariance matrix
psy_mat = noise2covar(n_mat); % psy_mat is the given noise covariance matrix
ncor_mat = n0_mat * chol(psy_mat) / sqrt(2); % generate noise data according to psy_mat 
psy_mat_new = noise2covar(ncor_mat); % recompute noise covariance matrix for checking
dpsy_mat = psy_mat - psy_mat_new; % difference matrix

% Construct signal+noise data for regenerated noise from Psy
sn_mat = zeros(nb_sample, nb_port); % initialize 
for k = 1:nb_port 
    sn_mat(:, k) = s_vec(k) + ncor_mat(:, k); % add noise to signals
end

% Compute direct sum-of-squares with regenerated noise
sos_dir = sos(sn_mat); % compute SoS
signal = mean(sos_dir); % evaluate signal
noise = std(sos_dir); % evaluate noise
snr_dir = signal / noise; % compute snr
fprintf('SNR from direct computation (regenerated noise) = %f\n', snr_dir);

disp(psy_mat);
disp(psy_mat_new);
disp(psy_mat-psy_mat_new);

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