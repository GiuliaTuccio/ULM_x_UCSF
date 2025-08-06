function data = PD_processing(PD,dataIQ)
%DOPPLER Summary of this function goes here
%   Detailed explanation goes here
[b,a] = butter(PD.B_order,PD.B_cutoff/(PD.fs/2),'high');

count=1;
data = [];

for i = 1:PD.n_e:size(dataIQ,3)-PD.n_e
    ens = dataIQ(:,:,i:i+PD.n_e-1);
    % subtract
    for k=1:size(ens,3)
    ens(:,:,k) = ens(:,:,k)-ens(:,:,1);
    end
    
    % Butter filter for tissue
    ens = filter(b,a,ens,[],3);
    
    % Discard first images
    ens = ens(:,:,PD.B_order+1:end);
    
    dataFlt = doppler_filter(ens,PD.B_svd);
    
    powerDoppler = mean(abs(dataFlt).^2,3);
    
    % bmode
    data(:,:,count) = 20*log10(abs(powerDoppler));
    count=count+1;
end
%%
if PD.plot_flag
    figure()
    for frame=1:size(data,3)
        imagesc(data(:,:,frame)-max(max(max(data))),[-80 0])
        colormap hot; drawnow; pause(0.5)
    end
end

