function [db_filter, db_background] = ULM_proc(ULM,dataIQ)
%ULM_PROC Summary of this function goes here
%   Detailed explanation goes here
n_ULM = 800;
s_idx = 1;
while s_idx <= ULM.nFrames
    end_idx = min(s_idx + n_ULM-1,ULM.nFrames);
    db_filter(:,:,s_idx:end_idx) = doppler_filter(dataIQ(:,:,s_idx:end_idx),ULM.k_value);
    if ULM.motion_compensation
        db_background(:,:,s_idx:end_idx) = inv_doppler(dataIQ(:,:,s_idx:end_idx),2);
    else
        db_background = zeros(size(dataIQ,1),size(dataIQ,2),s_idx:end_idx);
    end
    s_idx = end_idx + 1;
end

for f = 1:ULM.nFrames
    db_filter(:,:,f) = db_filter(:,:,f) - db_filter(:,:,1);
end

% %data_filter(~isfinite(data_filter))=0;
[b,a] = butter(ULM.B_order,ULM.B_cutoff/(ULM.fs/2),'high');
%Butter filter for tissue
db_filter = filter(b,a,db_filter,[],3);
end

