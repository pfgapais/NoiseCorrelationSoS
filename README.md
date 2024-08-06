# NoiseCorrelationSoS

CovPsySNRbasisSOS.m does not require additionnal data; the user can change the number of samples, the noise and signal properties or the number of coils.

CovPsySNRbasisSOS_WithDataInput requires a noise covariance matrix arrange in a (N x N) .mat file (N being the number of channels), and a coil sensitivity map arrange in a (n x N) .mat file (n being the number of voxels). The user can change the number of samples and the noise properties.
