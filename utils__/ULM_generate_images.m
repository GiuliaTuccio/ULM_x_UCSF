function [MatOut,MatOut_vel,MatOut_z] = ULM_generate_images(tracks_interp,resolution,sizeOut)
%ULM_GENERATE_IMAGES Summary of this function goes here
%   Detailed explanation goes here
aa=[0, 0, 0, 0]; % get the origin [ 1 1]
bb=[1, 1, 1, 1]*resolution; % get the slope
MatOut = zeros(sizeOut);
MatOut_vel = zeros(sizeOut);
MatOut_z = zeros(sizeOut);

TrackMatOut = cellfun(@(x) (x(:,[1 2 3 4])+aa).*bb, tracks_interp,'UniformOutput',0);
for i=1:size(TrackMatOut,2)

    % round pixel position into [z,x]
    pos_z_round=round(TrackMatOut{i}(:,1));
    pos_x_round=round(TrackMatOut{i}(:,2));

    % evaluate and smooth the velocity
    velnorm = smooth(vecnorm(TrackMatOut{i}(:,3:4),2,2),10);
    
    % if I want the axail direction, substitute velnorm with vel_z
    vel_z = velnorm.*sign(mean(TrackMatOut{i}(:,3)));

    % remove out of grid bubbles (ie. the grid is too small)
    outP = zeros(numel(pos_z_round),1);
    outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
    outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));
    
    ind=(sub2ind(sizeOut,pos_z_round(~outP),pos_x_round(~outP)));
    
    % tracks are counted only once per pixel
    ind=unique(ind);
    velnorm = velnorm(~outP);
    vel_z = vel_z(~outP);
    
    for j=1:size(ind,1)
        % position of the MBs
        MatOut(ind(j))=MatOut(ind(j))+1;
        % velocities
        MatOut_vel(ind(j)) = MatOut_vel(ind(j))+velnorm(j);
        % axis
        MatOut_z(ind(j)) = MatOut_z(ind(j))+vel_z(j);
    end
    
end
    

  
end

