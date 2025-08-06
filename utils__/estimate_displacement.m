function delta = estimate_displacement(ref_data,img2)
%ESTIMATE_DISPLACEMENT Summary of this function goes here
%   Detailed explanation goes here
% Get the size
[M, N] = size(img2);

% Compute FFTs
F1 = fft2(ref_data, M, N);
F2 = fft2(img2, M, N);

% Multiply F1 with complex conjugate of F2 (cross-correlation)
corr_freq = F1 .* conj(F2);

% Inverse FFT to get cross-correlation in spatial domain
corr_spatial = ifft2(corr_freq);

% Shift zero-frequency component to the center
corr_spatial_shifted = fftshift(real(corr_spatial));

% get center
[m,n] = find(max(max(corr_spatial_shifted))==corr_spatial_shifted);

%delta
delta = [m-M/2 n-N/2];

if isempty(delta)
    delta = [0 0];
end

