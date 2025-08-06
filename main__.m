clc; close all; clear all;
%% Ultra x UCSF in vivo data processing for ULM imaging

addpath(genpath('utils')) % add functions
%%
% 1. load bmf data

path = 'D:\UCSF\15 Jul\';
name_all = {'20250715T163228','20250715T163538','20250715T163943','20250715T164331'...
           '20250715T164728','20250715T165012','20250715T175355','20250715T175727','20250715T180026'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PD Info %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PD.B_order      = 10;         % Butter Order
PD.B_cutoff     = 50;         % Butter cutoff [Hz]
PD.B_svd        = 1;          % svd eigenvectors 
PD.n_e          = 200;        % ensemble dimension
PD.plot_flag    = false;      % plot flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ULM Info %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ULM.B_order = 10;             % Butter Order
ULM.B_cutoff = 60;            % [Hz]
ULM.k_value = 8;
ULM.motion_compensation = true;

t = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');

% Convert to string for folder name
folderName = [path filesep 'Results_' char(t)];


% Create the folder
mkdir(folderName);


for i=1:size(name_all,2)
name = [path filesep name_all{i}];
dirname = dir([name '\em*']);
% Extract numeric part from folder names
frameNames = {dirname.name};
folderNums = cellfun(@(x) sscanf(x, 'em_bmfData%d'), frameNames);
nFrames = size(dirname,1);
load([name '\postAcqParams.mat'])

[dataIQ, t] = read_data_ordered(folderNums,dirname,nFrames,Aux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([name filesep  'seqParams.mat']);
PD.fs = 1/Bmode.TprfEff;  % Sampling frequency [Hz]
PD_data = PD_processing(PD,dataIQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ULM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ULM.fs = PD.fs;
ULM.frameRate = nFrames/t(end);
ULM.nFrames = nFrames;

[ULM_data, ULM_background] = ULM_proc(ULM,dataIQ);

save([folderName filesep name_all{i}],'ULM_data',"ULM_background",'ULM','PD_data','PD')
end
%% -- 2. ULM processing---------
parpool('local', 4);
%%
folderName = 'D:\UCSF\15 Jul\Results_2025-07-21_12-35-50\talamus\';
ULM.nb = 100;
ULM.loc_method = 1; % imregionalmax
ULM.gap = 0;
ULM.n_ULM = 400;  % divide in blocks for efficiney
motion_compensation = [true false];
scaling_factor = [2 3];
linking = [2.5 5];
[a,b,c] = ndgrid(motion_compensation,scaling_factor,linking);
comb = [a(:), b(:), c(:)];

folder_db = dir([folderName filesep '2025*.mat']);
tracks_all = cell(1, size(folder_db, 1));  % Preallocate tracks_all

for c = 1:size(comb,1)
    ULM.motion_compensation = comb(c,1);
    ULM.scaling_factor = comb(c,2);
    ULM.linking = comb(c,3)*comb(c,2);

    ref_1 = load([folder_db(1).folder filesep folder_db(1).name]);
    ULM.PData = size(ref_1.ULM_data);
    ref_intra = imresize(abs(ref_1.ULM_background(:,:,1)),ULM.scaling_factor);
    ref_intra = 20*log10(abs(ref_intra));

parfor i = 1:size(folder_db, 1)
    % Load the data
    tmp = load([folder_db(i).folder filesep folder_db(i).name]); 
    frameRate = tmp.ULM.frameRate;   
    s_idx = 1;
    
    
    local_tracks = {};  % To hold tracks for each file
    
    while s_idx <= tmp.ULM.nFrames
        end_idx = min(s_idx + ULM.n_ULM - 1, tmp.ULM.nFrames);
        disp(['Start: ' num2str(s_idx) ', End: ' num2str(end_idx)])
        data = imresize(abs(tmp.ULM_data(:,:,s_idx:end_idx)),ULM.scaling_factor);
             
        % -inter-block
        current_frame = imresize(abs(tmp.ULM_background(:,:,s_idx)),ULM.scaling_factor);
        current_frame = 20*log10(abs(current_frame));
         
        % -Detection
        [~, points] = ULM_detect(data, ULM.nb, false);
        nFrames = max(points(:,3)) - min(points(:,3));
        
        if ULM.motion_compensation
            delta = estimate_displacement(ref_intra,current_frame);
            delta = [delta(1) -delta(2)];
            %data_motion = imresize(abs(tmp.ULM_background),scaling_factor);
            %data_motion = 20*log10(abs(data_motion));
            %points = motion_correction(data_motion,current_frame,points);
            points(:,1:2) = points(:,1:2) + repmat(delta, size(points,1),1);
            %points = compensate_MUST(ref_intra, current_frame, points);
        end
        
        % Tracking
        [~, tracks] = tracking(points, nFrames, ULM.linking, ULM.gap, 1, 1, frameRate, 2);
        
        % Accumulate tracks for this iteration
        local_tracks = cat(2, local_tracks, tracks);
        s_idx = end_idx + 1;
        
    end
    
    % Store the tracks for this file in tracks_all (safely inside parfor)
    tracks_all{i} = local_tracks;    
end

t = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
% Save the tracks after all iterations are completed
save([folderName filesep 'Tracks_' char(t)], 'tracks_all','ULM');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Generate the images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PData = [97 97];
lambda = 1.040540540540541e-01;
folderName = 'D:\UCSF\15 Jul\Results_2025-07-21_12-35-50\visual_cortex\';
xaxis = [-5:5];
zaxis = [0:10];
resolution = 3;
dir_tracks = dir([folderName filesep 'Tracks_*.mat']);
for i=1:size(dir_tracks,1)
    tracks = load([dir_tracks(i).folder filesep dir_tracks(i).name]);

    original_dim = PData*tracks.ULM.scaling_factor+[1 1];
    clear img img_v img_z
    for j=1:size(tracks.tracks_all,2)  
        %tracks = dir_tracks.tracks_all{j};    
        [img(:,:,j), img_v(:,:,j), img_z(:,:,j)] = ULM_generate_images(tracks.tracks_all{j},resolution,original_dim*resolution);
    end
    
    figure()
    imagesc(imgaussfilt(sum(img,3).^(1/3),0.8),[0 70])
    colormap hot

figure,
IntPower=1/3;SigmaGauss=0;
velColormap = cat(1,flip(flip(hot(128),1),2),hot(128)); % custom velocity colormap
velColormap = velColormap(5:end-5,:); % remove white parts
im=imagesc(sum(img_v,3).^IntPower.*sign(sum(img_z,3)));
im.CData=im.CData-sign(im.CData)/2;
colormap(gca,velColormap)
caxis([-1 1]*max(caxis)*.7) % add saturation in image
axis image; axis off

% Display MatOut Velocimetry
figure();
MatOut_vel = sum(img_v,3);
MatOut = sum(img,3);
clbsize = [180,50];
vmax_disp  = ceil(quantile(MatOut_vel(abs(MatOut_vel)>0),.98)/10)*10;
Mvel_rgb = MatOut_vel/vmax_disp; % normalization
Mvel_rgb(1:clbsize(1),1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))',1,clbsize(2)); % add velocity colorbar
Mvel_rgb = Mvel_rgb.^(1/1.5);Mvel_rgb(Mvel_rgb>1)=1;
Mvel_rgb = imgaussfilt(Mvel_rgb,.5);
Mvel_rgb = ind2rgb(round(Mvel_rgb*256),jet(256)); % convert ind into RGB

MatShadow = MatOut;MatShadow = MatShadow./max(MatShadow(:)*.3);MatShadow(MatShadow>1)=1;
MatShadow(1:clbsize(1),1:clbsize(2))=repmat(linspace(0,1,clbsize(2)),clbsize(1),1);
Mvel_rgb = Mvel_rgb.*(MatShadow.^IntPower);
Mvel_rgb = brighten(Mvel_rgb,.4);
BarWidth = round(1./(ULM.scaling_factor*lambda)); % 1 mm
Mvel_rgb(size(MatOut,1)-50+[0:3],60+[0:BarWidth],1:3)=1;
imshow(Mvel_rgb);axis on
title(['Velocity magnitude (0-' num2str(vmax_disp) 'mm/s)'])

end