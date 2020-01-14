%% calculate trasformation matrices
% @ D:\s_20190801_autoro_QB07_TIB3_01\fsems_prep_hum_rot_CWT1r_combi_rottrigger_0_180_10degstep_01.fid
clear all, close all, clc
tic
DATA = aedes_readfid(pwd); % import data
T2_trans = DATA.FTDATA(:,:,1:21:end);
T2_trans(:,94:end,:) = 3000 + randn(128,128-93,19)*2000;
T2_trans(T2_trans < 1.5e4) = 0;
ang = multiangle_coregister(T2_trans,0:10:180);

%% Calculate T2 maps

bg = double(mean(mean(DATA.FTDATA(3:35,3:79,1))));
T2 = zeros([128 128 19]);
TE = DATA.PROCPAR.pw_mlev(1:6)*1e-3;
TE = TE(:) - TE(1);
[m,n] = size(T2(:,:,1));
thres = 15e3;
for NN = 1:19
    ind = ((NN-1)*21+1);
    T2_raw{NN}.FTDATA = DATA.FTDATA(:,:,ind:ind+5);
    T2_raw{NN}.PROCPAR = DATA.PROCPAR;
    T2_raw{NN}.PROCPAR.te = TE;
    T2_raw{NN}.HDR = DATA.HDR;
    
    h=waitbar(0,['calculating ', int2str(NN), '/19']);
    drawnow;
    for ii=1:m
        for jj=1:n
            s=double(squeeze(T2_raw{NN}.FTDATA(ii,jj,:)));
            if(sum(s) ~= 0)
                T2dat=calculate_t2(TE,s,2,0,bg);
                if (T2dat(2) > thres)
                    T2dat = [0 0 0];
                end
            else
                T2dat=[0 0 0];
            end
            M0(ii,jj)=T2dat(1);
            T2(ii,jj,NN)=T2dat(2);
            
            Const(ii,jj)=T2dat(3);
            %Err(ii,jj)=T2dat(4);
        end
        waitbar(ii/m,h);
    end
    close(h);
    
end
%% anisotropy

multiangle_transformix(T2,0:10:180);
T2_reged = aedes_read_nifti('Registered_maps.nii');
for ii = 1:19
    T2_reged.FTDATA(:,:,ii) = T2_reged.FTDATA(:,:,ii).*ROI_mask.voxels{1};
end
aedes_write_nifti(T2_reged,'Registered_T2_maps.nii')
%%
T2_reged = aedes_read_nifti('Registered_T2_maps.nii');

R2 = T2_reged.FTDATA.^-1;
R2(isnan(R2)) = 0;
R2(isinf(R2)) = 0;
ma = max(R2,[],3);
mi = min(R2,[],3);
anisotropymap = (ma-mi)./(ma+mi);
anisotropymap(isnan(anisotropymap)) = 0;
anisotropymap(anisotropymap > 1) = 1;
anisotropymap(anisotropymap < 0) = 0;
anisocrop = anisotropymap([79:191],[52:178],:);
imagesc(anisocrop)
colorbar
colormap jet

%% R2 vs angle profile for zone avgs

temp = load('ROI_zones.roi','-mat');
ROI_zones = temp.ROI;


for lay_m = 1:19
    for ii = 1:3
        [row_m{ii},col_m{ii}] = find(ROI_zones(ii).voxels{1});
        lin_ind = sub2ind([128 128 19],row_m{ii}, col_m{ii}, lay_m*ones(length(col_m{ii}),1));
        zone_avg(lay_m,ii) = mean(T2_reged.FTDATA(lin_ind));
    end
end

% Gauss-Newton iteration
th = [0.02 0.02]';
er = 1;
ite = 1;
step = 1;
th(3) = -2;
while ((norm(er) > 1e-6) && (ite < 10000)) % jaakoppi
    J(:,1) = 1*ones(length(ang),1) ;
    J(:,2) = (3*cosd(ang+th(3)).^2 -1).^2;
    J(:,3) = 2*th(2)*((cosd(ang+th(3))).^2 - 1).*2.*cosd(ang+th(3)).*sind(ang+th(3));
    %J(:,1) = 0*ones(length(ang),1) ;
    r = zone_avg(:,3).^-1 - (th(1) + th(2)*(3*cosd(ang + th(3)).^2-1).^2);
    thnew = th + step*inv(J'*J)*J'*r;
    er = abs(thnew-th);
    th = thnew;
    ite = ite+1;
    teetta(:,ite) = th;
end
th
ite
%th(3) = -5;

figure
hold on
%plot(ang,zone_avg,'*-')
plot(ang,zone_avg(:,1).^-1,'o-','Color',[0.4,0.7,0.9],'linewidth',1.5)
plot(ang,zone_avg(:,2).^-1,'o-','Color',[0 0.3 1],'linewidth',1.5)
plot(ang,zone_avg(:,3).^-1,'o-','Color',[0 0 0],'linewidth',1.5)
plot(0:180, th(1) + th(2)*(3*cosd([0:180]+th(3)).^2-1).^2,'r--','linewidth',1.5)
legend('SZ','TZ','RZ','R_{2,i} + R_{2,a}(3cos^2(\theta + \phi) - 1)^2')
xlabel('Kulma ( ^{\circ} )')
ylabel('R_2 (Hz)')

%% T2 vs depth

for kulma = 1:19
    plotti = new_roi2profile(T2_reged.FTDATA(:,:,kulma),ROI_zones,14);
    depthprofiles(:,kulma) = equisample(plotti,100);
end
[X,Y] = meshgrid(1:100,ang);
figure
hold on
surf(Y,X,depthprofiles')
xlabel('Kulma (deg)')
ylabel('Syvyys (%)')
zlabel('T_2 (ms)')

%% Phasemaps

load ROIs_regd.mat
T2_reged = aedes_read_nifti('Registered_T2_maps.nii');

phasemap = zeros(size(ROI_mask.voxels{1}));
h=waitbar(0,['calculating ']);
disced = 0;
for ii = 1:size(T2_reged.FTDATA,1)
    for jj = 1:size(T2_reged.FTDATA,2)
        waitbar(ii/128,h,[num2str(ii), '/' , num2str(128)]);
        
        if (ROI_mask.voxels{1}(ii,jj))
            th = [0.02 0.02]';
            th(3) = -2;
            er = 1;
            ite = 1;
            step = 1;
            while ((norm(er) > 1e-6) && (ite < 10000)) 
                J(:,1) = 1*ones(length(ang),1) ; % jaakoppi
                J(:,2) = (3*cosd(ang+th(3)).^2 -1).^2;
                J(:,3) = 2*th(2)*((cosd(ang+th(3))).^2 - 1).*2.*cosd(ang+th(3)).*sind(ang+th(3));
                %J(:,1) = 0*ones(length(ang),1) ;
                r = squeeze(T2_reged.FTDATA(ii,jj,:)).^-1 - (th(1) + th(2)*(3*cosd(ang + th(3)).^2-1).^2);
                thnew = th + step*inv(J'*J)*J'*r;
                er = abs(thnew-th);
                th = thnew;
                ite = ite+1;
                if (norm(th(3)) > 180 || isnan(norm(th)))
                    disced = disced+1;
                    th = [0 0 0]';
                    continue;
                end
            end
            phasemap(ii,jj,1) = th(1);
            phasemap(ii,jj,2) = th(2);
            phasemap(ii,jj,3) = th(3);
        end
    end
end
disp(disced)
close(h);

%%
DATA = aedes_readfid(pwd);
bg = double(mean(mean(DATA.FTDATA(3:35,3:79,1))));
DATA = flip(flip(DATA.FTDATA,1),2);
stecho = DATA(:,:,1:6);
orientationsT2 = DATA(:,:,1:21:end);

kaura = double(stecho(30,62,:));
kaura = kaura(:);
TE = [0 8 16 32 64 128]';
T2dat=calculate_t2(TE,kaura,3,0,bg);

M = T2dat(1).*exp(-[0:128]'./T2dat(2)) + T2dat(3);
figure
plot(TE,kaura,'ro','linewidth',1.5);
hold on
plot([0:128]',M,'linewidth',1.5);
xlabel('Aika (ms)')
ylabel('Magnitudi')
legend('Havainnot','Sovitettu käyrä M = M_0 e^{-t/T_2} + C')














