function doppler_data = doppler_filter(data,n_svd)
%DOPPLER_FILTER Summary of this function goes here
%   Detailed explanation goes here
data(isnan(data)) = 0;
[nz,nx,nt] = size(data);

% rearrange
m = zeros(nx*nz,nt);
i=0;
for ix = 1:nx
    for iz = 1:nz
        i=i+1;
        m(i,:) = squeeze(data(iz,ix,:)-mean(data(iz,ix,:)));
    end
end

[U,S,V] = svd(m,'econ');

Sf = S;
Sf(1:n_svd,1:n_svd) = 0;
mf = U*Sf*V';

i=0;
doppler_data = zeros(nz,nx,nt);

for ix=1:nx
    for iz=1:nz
        i=i+1;
        doppler_data(iz,ix,:) = squeeze(mf(i,:));
    end
end