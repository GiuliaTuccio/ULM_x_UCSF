function [dataIQ,t] = read_data_ordered(folderNums,dirname,nFrames,Aux)
%READ_DATA_ORDERED Summary of this function goes here
%   Detailed explanation goes here
% Sort based on numeric value
[~, sortIdx] = sort(folderNums);
sortedFolders = dirname(sortIdx);

t_0 = Aux.frameRate{1};

disp([sortedFolders(1).folder])

for i = 1:nFrames
    strs = strsplit(sortedFolders(i).name,'em_bmfData');
    %load([sortedFolders(i).folder filesep sortedFolders(i).name filesep 'em' strs{2} '.mat'])
    load([sortedFolders(i).folder filesep sortedFolders(i).name ])
    
    t_now = Aux.frameRate{i};
    t(i) = seconds(t_now-t_0);
    %imagesc(bimg)
    %drawnow; pause(0.1)
    %data(:,:,i) = double(bimg(:,:));
    dataIQ(:,:,i) = double(bmfData(:,:));
end
end

