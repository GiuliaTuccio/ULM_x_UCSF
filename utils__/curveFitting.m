function [z_c, x_c, sigma]=curveFitting(I_ROI,vectwhmz,vectwhmx)
    [mesh_x, mesh_z]=meshgrid(vectwhmx,vectwhmz);
    meshIn=cat(3,mesh_x,mesh_z);

    sigma_z = 1;
    sigma_x = 1;
    myGaussFun = @(x_pos, mesh_pos)(exp(-(mesh_pos(:,:,1)-x_pos(1)).^2./(2*sigma_z^2)-(mesh_pos(:,:,2)-x_pos(2)).^2./(2*sigma_x^2)));
    OPTIONS=optimoptions('lsqcurvefit','StepTolerance',.01,'MaxIterations',10,'Display','off');

    % fit
    x_out=lsqcurvefit(myGaussFun, [0,0], meshIn, double(I_ROI./(max(I_ROI(:)))), [], [], OPTIONS);
    z_c = x_out(2);
    x_c = x_out(1);
    sigma = ComputeSigmaScat(I_ROI,z_c,x_c);
end
function sigma = ComputeSigmaScat(Iin,Zc,Xc)
%% This function will calculate the Gaussian width of the presupposed peak in the intensity, which we set as an estimate of the width of the microbubble
    [Nx,Nz] = size(Iin);
    Isub = Iin - mean(Iin(:));
    [px,pz] = meshgrid(1:Nx,1:Nz);
    zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;
    xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
    r2 = zoffset.*zoffset + xoffset.*xoffset;
    sigma = sqrt(sum(sum(Isub.*r2'))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
end