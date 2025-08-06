function MatTracking = localization(data, mask_idx_x,mask_idx_z,mask_frames, fwhmz, fwhmx,d)
        %LOCALIZATION Summary of this function goes here
        %   Detailed explanation goes here
        
        if nargin<5 || isempty(fwhmz)
            fwhmz=3;
        end
        if nargin<6 || isempty(fwhmx)
            fwhmx=3;
        end

        [width, height, ~] = size(data);
        
        
        % initialize the results
        avarageXc = nan(size(mask_idx_x,1),1);
        avarageZc = nan(size(mask_idx_z,1),1);
        for nb=1:size(mask_idx_z,1)
            vectwhmz = -1*floor(fwhmz/2):floor(fwhmz/2);
            vectwhmx = -1*floor(fwhmx/2):floor(fwhmx/2);
           % select the region of interest
        %        if mask_idx_z(nb)+vectwhmz(end) > width || mask_idx_x(nb)+vectwhmx(end) > height
        %            disp(nb)
        %            I_ROI = data(mask_idx_z(nb)+vectwhmz(1):width, mask_idx_x(nb)+vectwhmx(1):height,mask_frames(nb));
        %            if(nnz(d(mask_idx_z(nb)+vectwhmz(1):width, mask_idx_x(nb)+vectwhmx(1):height,mask_frames(nb))))>1
        %             continue
        %            end 
        %            vectwhmz=mask_idx_z(nb)+vectwhmz(1):width;
        %            vectwhmx=mask_idx_x(nb)+vectwhmx(1):height;
        %        elseif mask_idx_z(nb)+vectwhmz(1) < 1 || mask_idx_x(nb)+vectwhmx(1) < 1
        %             I_ROI = data(1:mask_idx_z(nb)+vectwhmz(end), 1:mask_idx_x(nb)+vectwhmx(end),mask_frames(nb));

        %            vectwhmz=1:mask_idx_z(nb)+vectwhmz(end);
        %            vectwhmx=1:mask_idx_x(nb)+vectwhmx(end);
        %        else
        %            I_ROI = data(mask_idx_z(nb)+vectwhmz, mask_idx_x(nb)+vectwhmx,mask_frames(nb));
        %            if(nnz(d(mask_idx_z(nb)+vectwhmz, mask_idx_x(nb)+vectwhmx,mask_frames(nb))))>1
        %            continued
        %             end
        %        end 
        %    disp(nb)
            start_z = mask_idx_z(nb)+vectwhmz(1);
            start_x = mask_idx_x(nb)+vectwhmx(1);
            end_z = mask_idx_z(nb)+vectwhmz(end);
            end_x = mask_idx_x(nb)+vectwhmx(end);
            if end_z > width 
                end_z = width;
            elseif start_z   < 1
                start_z = 1;
            end 
            if end_x > height 
                end_x = height;
            elseif start_x   < 1
                start_x = 1;
            end 
            I_ROI = data(start_z:end_z, start_x:end_x,mask_frames(nb));
            if(nnz(d(start_z:end_z, start_x:end_x,mask_frames(nb))))>1
            continue
           end
           % localization method through Gaussian fit
           [z_c, x_c, sigma]=curveFitting(I_ROI,start_z-mask_idx_z(nb):end_z-mask_idx_z(nb),start_x-mask_idx_x(nb):end_x-mask_idx_x(nb));
        
           % save the final result
           avarageZc(nb)=z_c+mask_idx_z(nb);
           avarageXc(nb)=x_c+mask_idx_x(nb);
        
           % Additional safeguards
           % sigma evaluates the size of the microbubble. If it appears to be too large, the microbubble can be removed (optional)
           if or(sigma<0,sigma>25)
           %         averageZc(iscat)=nan;
           %         averageXc(iscat)=nan;
           %         continue
           end
            
           % If the final axial/lateral shift is higher that the fwhmz,
           % localization has diverged and the microbubble is ignored.
           if or(abs(z_c)>fwhmz/2,abs(x_c)>fwhmx/2)
               avarageZc(nb)=nan;
               avarageXc(nb)=nan;
               continue
           end
        end
        
        % check for null
        keepIdx = ~isnan(avarageXc);
        
        % save the results
        MatTracking = zeros(nnz(keepIdx),3);
        MatTracking(:,1)=avarageZc(keepIdx);
        MatTracking(:,2)=avarageXc(keepIdx);
        MatTracking(:,3)=mask_frames(keepIdx);

        clear avarageZc avarageXc mask_frames keepIdx
        end

