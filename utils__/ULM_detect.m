function [MaskFinal, points] = ULM_detect(data, noMB, loc)
%mask_idx_z,mask_idx_x, mask_frame, MaskFinal
% LOCALIZATION function that performs the detection of the microbubbles,
% given an intensity matrix as input
%  Input:
%        data: matrix/tensor containing the intensity for each pixel. Input
%        must be a 3d matrix (z, x, nFrames)
%        method: 
%              1: "imregionalmax" MATLAB as default method. 
%              2:  kernel method, which is the same as imregional max, but considers
%                  a kernel in which to find the local maxima. If "kernel" a kernel
%                  size must be provided. kernel must be even number!!
%        noMB: number of microbubbles
%
%  Output:
%        mask_idx_z: vector containing the index of the z coordinates of
%        the detections
%        mask_idx_x: vector containing the index of the x coordinates of
%        the detections
%        mask_frame: vector containing the index of the frame of
%        the detections

   
    if nargin<2 || isempty(noMB)
        disp('Provide an approximative number of MB')
    end
    if size(data,2) < 1 || isempty(data)
        disp('Input data must have at least 1 frame')
    end
    if nargin<3
        loc = false;
    end

    

    [height,width,numberOfFrames] = size(data);
    % pre-allocate memory
    mask2D = zeros(height,width,numberOfFrames);

    % save the highest pixel value in a mask
    mat2D  = permute(data, [1,3,2]);
    mat2D  = reshape(mat2D,height*numberOfFrames,width);
    mask   = imregionalmax(mat2D); clear mat2D
    mask2D = reshape(mask, height,numberOfFrames,width); clear mask
    mask2D = permute(mask2D,[1,3,2]); % so that we restore (z,x,t) table
    
    
    % keep just the strongest signal
    IntensityMatrix_new = data.*mask2D;
    [tempMat,~] = sort(reshape(IntensityMatrix_new,[],size(IntensityMatrix_new,3)),1,'descend');
    IntensityFinal = IntensityMatrix_new - ones(size(IntensityMatrix_new)).*reshape(tempMat(noMB,:),[1 1 numberOfFrames]);
    
    MaskFinal = (mask2D.*IntensityFinal)>0;
    clear mask2D IntensityMatrix_new tempMat
    MaskFinal(isnan(MaskFinal))=0;
    MaskFinal = (MaskFinal>0).*IntensityFinal;

    % hold the information
    [mask_idx_z,mask_idx_x,mask_frame]= ind2sub([height,width,numberOfFrames],find(MaskFinal));
    

    if loc 
        points = localization(data,mask_idx_x,mask_idx_z,mask_frame,3,3,MaskFinal);
    else
        points = cat(2, mask_idx_z, mask_idx_x, mask_frame);
    end
    
    
    clear IntensityFinal 
end

   
