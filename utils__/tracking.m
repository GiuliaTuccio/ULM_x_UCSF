function [MatTracking,tracks_out] = tracking(MatTracking, nFrames, max_linking_distance, max_gap, dx, dz, frameRate, track_length,interp)
    
    % matTracking : 
    %               1. x axis
    %               2. z axis
    %               3. frame_id
    
    % method: tracking method: 0: simple tracking based on the bipartite
    %                              matching
    %                          1: kalman filtering
    %

    if nargin<8
        track_length = 2;
    end
    if nargin<9
        interp = true;
    end


    % retrieve information of the points
    frame_idx = arrayfun(@(i) find(MatTracking(:,3)==i), [1:nFrames],'UniformOutput',false);
    Points = arrayfun(@(i) [MatTracking(frame_idx{i},1),MatTracking(frame_idx{i},2)],[1:nFrames],'UniformOutput',false);
    
    tracks_raw = {};
    
   
%     figure()
%     for f=1:nFrames
%         drawnow; 
%         hold on
%         scatter(Points{f}(:,2),Points{f}(:,1),'x')
%     end
%     
    % simple_tracks is a cell per track found  Each track is made of a |n_frames x 1| integer array, containing the index of the particle belonging to that track in the corresponding frame. NaN values report that for this track at this frame, a particle
    % could not be found (gap).
    % adjency track return also a cell array with one cell per track, but the indices in each track are the global indices of the concatenated points array
    %tracks_raw = cell(1,5000);
    
        
    [simple_tracks, adjency_tracks]=simpletracker( Points, ...
        'MaxLinkingDistance', max_linking_distance, ...
        'MaxGapClosing', max_gap);
    
    % initialize number of tracks
    n_tracks=numel(simple_tracks);
    colors=hsv(n_tracks);
    % retrieve info about the points
    all_points=vertcat(Points{:});
    %used to count the tracks
    count=1;
    
    %figure()
    %hold on
    for t=1:n_tracks
        % retrieve each track
        track=adjency_tracks{t};


        % retreive the id frame
        idFrame=MatTracking(track,3);
        % retrieve the points id for th track
        track_points=cat(2,all_points(track,:),idFrame);
        % save if the length is shorter than ...
        if length(track_points(:,1))>2
            tracks_raw{count}=track_points.*[dz dx 1];
            % plot(tracks_raw{count}(:,2),tracks_raw{count}(:,1),'Color',colors(t,:))
            count= count+1;
        end
    end

    % interpolate
    interp_factor = 1/2;
    tracks_out = {};
    
    if interp
        
        for t=1:size(tracks_raw,2)
            if size(tracks_raw{t},1) > track_length
            zi = tracks_raw{t}(:,2);
            xi = tracks_raw{t}(:,1);

            TimeABS = (0:(length(zi)-1))*1/frameRate;

            z_interp = interp1(1:length(zi), smooth(zi, 50),1:interp_factor:length(zi));
            x_interp = interp1(1:length(xi), smooth(xi, 50),1:interp_factor:length(xi));
            t_interp=interp1(1:length(TimeABS),TimeABS,1:interp_factor:length(TimeABS));

            % Velocity
            vzu = diff(z_interp)./diff(t_interp); vzu = [vzu(1), vzu];
            vxu = diff(x_interp)./diff(t_interp); vxu = [vxu(1), vxu];

            tracks_out{t} = cat(2,x_interp',z_interp', vxu', vzu',t_interp');
            end
        end
    else
        tracks_out = tracks_raw;
    end

     tracks_out = tracks_out(~cellfun('isempty',tracks_out));

    clear track_points
   
     