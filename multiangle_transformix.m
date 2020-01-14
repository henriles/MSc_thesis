function f = multiangle_transformix(Stack,init_angle,folder)
% (c) Henri Leskinen, University of Eastern Finland, 2019
%
%Rotates image files (matrices) using transform transformation parameter-
% files (Transpara[init_angle].txt) genereted by multiangle_coregister.m
% -function. 
% Takes the images input as a 3D-matrix (Stacked images)
% Uses the first image as reference 


if (nargin < 3)
    folder = [pwd, '\Transform'];
end

for ii = 1:size(Stack,3)
    
    data = imrotate(Stack(:,:,ii),init_angle(ii),'crop','bilinear');
    aedes_write_nifti(data,[folder, '\Rotationdata_', int2str(init_angle(ii)), 'deg.nii']);
    
    if (ii == 1)
        % dont transform the first one
        tempdata = aedes_read_nifti([folder, ['\Rotationdata_', int2str(init_angle(1)), 'deg.nii']]);
        regstack = tempdata;
    else
        
        transformixstring = ['C:\Users\henriles\Documents\elastix-4.9.0-win64\transformix ',...
            '-in ', folder, '\Rotationdata_', int2str(init_angle(ii)), 'deg.nii ',....
            '-out ', folder, ' ',...
            '-tp ', folder, '\Transpara', int2str(init_angle(ii)), '.txt'];
        system(transformixstring);
        
        delete([folder, '\Rotationdata_', int2str(init_angle(ii)), 'deg.nii']);
        tempdata = aedes_read_nifti([folder, '\result.nii']);
        regstack.FTDATA(:,:,ii) = tempdata.FTDATA(:,:);
        
    end
    aedes_write_nifti(regstack,[folder, '\..\Registered_maps.nii']);
end
