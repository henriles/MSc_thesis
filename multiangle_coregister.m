function ang_co = multiangle_coregister(Stack,angle_guess,folder)
% (c) Henri Leskinen, University of Eastern Finland, 2019
%
% Coregisters images (matrices) to the first one in the stack. 
% Uses initial guess for the images angles offset to index files and 
% increase the propability of succesful registeration.
% Requires Elastix (http://elastix.isi.uu.nl) to be installed in order to
% work. 
% Uses image input as 3D matrix (Stacked images)

if (nargin < 3)
    folder = [pwd, '\Transform'];
end
if ~exist(folder, 'dir')
    mkdir(folder)
end
h=waitbar(0,['calculating ']);
drawnow;
for ii = 1:size(Stack,3)
    
    data = imrotate(Stack(:,:,ii),angle_guess(ii),'crop','bilinear');
    aedes_write_nifti(data,[folder, '\Rotationdata_', int2str(angle_guess(ii)), 'deg.nii']);
    
    if (ii == 1)
        tempdata = aedes_read_nifti([folder, '\Rotationdata_', int2str(angle_guess(1)), 'deg.nii']);
        regstack = tempdata;
    else
        elastringi = ['C:\Users\henriles\Documents\elastix-4.9.0-win64\elastix ',...
            '-f ', folder, '\Rotationdata_', int2str(angle_guess(1)), 'deg.nii ',...
            '-m ', folder, '\Rotationdata_', int2str(angle_guess(ii)), 'deg.nii ',...
            '-out ', folder, ' '...
            '-p ', folder, '\Parameters_Rigid.txt'];
        system(elastringi);
        tempdata = aedes_read_nifti([folder, '\result.0.nii']);
        regstack.FTDATA(:,:,ii) = tempdata.FTDATA;
        try
            delete([folder, '\Transpara', int2str(angle_guess(ii)), '.txt']);
        catch
            disp('Keijo');
        end
        system(['ren ', folder, '\TransformParameters.0.txt ',...
            'Transpara', int2str(angle_guess(ii)), '.txt']);
        
        nimi = [folder, '\Transpara', int2str(angle_guess(ii)), '.txt'];
        temp0 = importdata(nimi);
        fid=fopen(nimi, 'r');
        temp1 = textscan(fid, '%s');
        fclose(fid);
        ang_co(ii,1) = angle_guess(ii) - str2double(temp1{1}{6})*180/pi;
        waitbar(ii/size(Stack,3),h);
        delete([folder, '\Rotationdata_', int2str(angle_guess(ii)), 'deg.nii']);
    end
end
close(h);
aedes_write_nifti(regstack,[folder, '\Regstack.nii']);

end

