% @ D:\s_20190801_autoro_QB07_TIB3_01
% Calculate T2 maps from fid
clear all, close all, clc
DATA = aedes_readfid('fsems_prep_hum_rot_CWT1r_combi_rottrigger_0_180_10degstep_01.fid');
load 'fsems_prep_hum_rot_CWT1r_combi_rottrigger_0_180_10degstep_01.fid\T2_raw_and_ROIs.mat';

bg = double(mean(mean(DATA.FTDATA(3:35,3:79,1))));
T2 = zeros([128 128 19]);
TE = DATA.PROCPAR.pw_mlev(1:6)*1e-3;
TE = TE(:) - TE(1);
[m,n] = size(T2(:,:,1));
thres = 300;
for NN = 1:19
    ind = ((NN-1)*21+1);
    T2_raw{NN}.FTDATA = ROI.voxels{NN}.*DATA.FTDATA(:,:,ind:ind+5);
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
AddInfo.save = false;
save('T2_rot.mat','T2')

% ROI = [];
% for nn = 1:19
%     ROI.voxels{nn} = eval(['ROI_', int2str(nn), '.voxels{1}'])
% end

%%
%close all
%% calc angles and plot zone averages
clear all
load T2_rot.mat
load ROI_zones.mat

for ii = 1:4
    [row(:,ii),col(:,ii),lay(:,ii)] = ind2sub(size(ROI_zones(ii+3).voxels{1}),find(ROI_zones(ii+3).voxels{1}));
end
ang = atand((row(:,3) - row(:,4)) ./ (col(:,3)-col(:,4)));
ang(ang<0) = ang(ang<0) + 180;
%ang = ang-5;

for lay_m = 1:19
    for ii = 1:3
        [row_m{ii},col_m{ii}] = find(ROI_zones(ii).voxels{1}(:,:,lay_m));
        lin_ind = sub2ind([128 128 19],row_m{ii}, col_m{ii}, lay_m*ones(length(col_m{ii}),1));
        zone_avg(lay_m,ii) = mean(T2(lin_ind));
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
th(3) = -5;

figure
hold on
%plot(ang,zone_avg,'*-')
plot(ang,zone_avg.^-1,'*-')
plot(0:180, th(1) + th(2)*(3*cosd([0:180]+th(3)).^2-1).^2)
legend('SZ','TZ','RZ','C_1 + C_2(3cos^2(\theta + C_3)-1)^2')
xlabel('Angle (deg)')
ylabel('R2')

%% plot depth-profile vs angle
clear ROI

figure
hold on
ROI = ROI_zones;
for kulma = 1:19
    for N =1:7
        ROI(N).voxels{1} = [];
        ROI(N).voxels{1} = ROI_zones(N).voxels{1}(:,:,kulma);
    end
    All_rois{kulma} = ROI;
    plotti = new_roi2profile(T2(:,:,kulma),All_rois{kulma},14);
    A(:,kulma) = equisample(plotti,100);
end
[X,Y] = meshgrid(1:100,ang);
surf(Y,X,A')
%%
%% T1: construct a stacked version
clear all, close all, clc

for ii = 1:19
    filename = ['rot_', int2str((ii-1)*10), '\fsems_prep_hum_rot_IR_T1FSE_',...
        int2str((ii-1)*10), '_01.fid\t1_ir_map_2par.mat'];
    load(filename)
    DATA(:,:,ii) = Data;
end

%% T1
clear all
load ROI_zones_T1.mat
load T1_maps_stacked.mat
T1 = DATA.FTDATA;

for ii = 1:4
    [row(:,ii),col(:,ii),lay(:,ii)] = ind2sub(size(ROI_zones_T1(ii+3).voxels{1}),find(ROI_zones_T1(ii+3).voxels{1}));
end
ang = atand((row(:,3) - row(:,4)) ./ (col(:,3)-col(:,4)));
ang(ang<0) = ang(ang<0) + 180;
%ang = ang-5;

for lay_m = 1:19
    for ii = 1:3
        [row_m{ii},col_m{ii}] = find(ROI_zones_T1(ii).voxels{1}(:,:,lay_m));
        lin_ind = sub2ind([128 128 19],row_m{ii}, col_m{ii}, lay_m*ones(length(col_m{ii}),1));
        zone_avg(lay_m,ii) = mean(T1(lin_ind));
    end
end

% Gauss-Newton iteration
th = [0.02 0.02]';
er = 1;
ite = 1;
step = 1e-1;
while ((norm(er) > 1e-6) && (ite < 10000)) % jaakoppi
    J(:,1) = 1 + th(2)*(3*cosd(ang).^2 -1).^2;
    J(:,2) = th(1) + (3*cosd(ang).^2 -1).^2;
    r = zone_avg(:,3).^-1 - (th(1) + th(2)*(3*cosd(ang).^2-1).^2);
    thnew = th + step*inv(J'*J)*J'*r;
    er = abs(thnew-th);
    th = thnew;
    ite = ite+1;
    teetta(:,ite) = th;
end
th
ite

figure
hold on
%plot(ang,zone_avg,'*-')
plot(ang,zone_avg.^-1,'*-')
%plot(0:180, th(1) + th(2)*(3*cosd(0:180).^2-1).^2)
legend('SZ','TZ','RZ','C_1 + C_2(3cos^2(\theta)-1)^2')
xlabel('Angle (deg)')
ylabel('R1')


%% T1: plot depth-profile vs angle
clear ROI

figure
hold on
ROI = ROI_zones_T1;
for kulma = 1:19
    for N =1:7
        ROI(N).voxels{1} = [];
        ROI(N).voxels{1} = ROI_zones_T1(N).voxels{1}(:,:,kulma);
    end
    All_rois{kulma} = ROI;
    plotti = new_roi2profile(T1(:,:,kulma),All_rois{kulma},14);
    A(:,kulma) = equisample(plotti,100);
end
figure
[X,Y] = meshgrid(1:100,ang);
surf(Y,X,A')



%%
for ii = 2:18
    diffi(ii) = ang_co(ii+1)-ang_co(ii)
end

%%
clear all, close all, clc
tic
DATA = aedes_readfid(pwd);
T2_trans = DATA.FTDATA(:,:,1:21:end);
T2_trans(:,94:end,:) = 3000 + randn(128,128-93,19)*2000;
T2_trans(T2_trans < 1.5e4) = 0;
ang = multiangle_coregister(T2_trans,0:10:180);



multiangle_transformix(T2_trans,0:10:180);
toc




