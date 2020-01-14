%% read and stack T1 fids
clear all, aedes_killfigs, clc

for ii = 1:19
    filename = ['rot_', int2str((ii-1)*10), '\fsems_prep_hum_rot_IR_T1FSE_',...
        int2str((ii-1)*10), '_01.fid'];
    Data{ii} = aedes_readfid(filename);
    ind = (ii-1)*6 + 1;
    if ii == 1
        DATA = Data{ii};
    end
    DATA.FTDATA(:,:,ind:ind+5) = Data{ii}.FTDATA;
end

Data_th = DATA.FTDATA(:,:,1:6:end-6);
Data_th(Data_th < 10e3) = 0;
Data_th(10:end,100:end,:) = 0;
Data_th = flip(Data_th,3);

ang = multiangle_coregister(Data_th,0:-10:-170);
%% Calcucate T1 maps

T1 = zeros([128 128 18]);
ti = DATA.PROCPAR.ir_delay*1e3;
[m,n] = size(T1(:,:,1));
thres = 15e3;
for NN = 1:18
    
    h=waitbar(0,['calculating ', int2str(NN), '/19']);
    drawnow;
    for ii=1:m
        for jj=1:n
            s=double(squeeze(Data{19 - NN}.FTDATA(ii,jj,:)));
            if(sum(s) ~= 0)
                T1dat=calculate_t1ir(ti,s,2);
%                 if (T1dat(1) > thres)
%                     T1dat = [0 0];
%                 end
            else
                T1dat=[0 0];
            end
            M0(ii,jj,NN)=T1dat(1);
            T1(ii,jj,NN)=T1dat(2);
        end
        waitbar(ii/m,h);
    end
    close(h);
end

%% Use transformix to orient T1 maps. then apply mask (calcroi)

T1(T1 > 3000) = 0;
multiangle_transformix(T1,0:-10:-170);
T1_reged = aedes_read_nifti('Registered_maps.nii');
for ii = 1:18
    T1_reged.FTDATA(:,:,ii) = T1_reged.FTDATA(:,:,ii).*ROI_calc.voxels{1};
end
aedes_write_nifti(T1_reged,'Registered_maps.nii')
aedes(T1_reged)

%% R1 anisotropymap

T1_reged = aedes_read_nifti('Registered_maps.nii');

R1 = T1_reged.FTDATA.^-1;
R1(isnan(R1)) = 0;
R1(isinf(R1)) = 0;
ma = max(R1,[],3);
mi = min(R1,[],3);
anisotropymap = (ma-mi)./(ma+mi);
anisocrop = anisotropymap([43:116],[27:87],:);
anisotropymap(isnan(anisotropymap)) = 0;
anisotropymap(anisotropymap > 1) = 1;
anisotropymap(anisotropymap < 0) = 0;
aedes(anisotropymap)

%% R1 vs angle profile for zone avgs
for lay_m = 1:size(T1_reged.FTDATA,3)
    for ii = 1:3
        [row_m{ii},col_m{ii}] = find(ROI_zones(ii).voxels{1});
        lin_ind = sub2ind([128 128 19],row_m{ii}, col_m{ii}, lay_m*ones(length(col_m{ii}),1));
        zone_avg(lay_m,ii) = mean(T1_reged.FTDATA(lin_ind));
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
%ang = -ang + 12.5288;
zone_avg = zone_avg/1000;
figure
hold on
axis([0 180 0 8e-1])
%plot(ang,zone_avg,'*-')
plot(ang,zone_avg(:,1).^-1,'o-','Color',[0.4,0.7,0.9],'linewidth',1.5)
plot(ang,zone_avg(:,2).^-1,'o-','Color',[0 0.3 1],'linewidth',1.5)
plot(ang,zone_avg(:,3).^-1,'o-','Color',[0 0 0],'linewidth',1.5)
%plot(0:180, th(1) + th(2)*(3*cosd([0:180]+th(3)).^2-1).^2,'r')
legend('SZ','TZ','RZ');%,'C_1 + C_2(3cos^2(\theta + C_3) - 1)^2')
xlabel('Kulma ( ^{\circ} )')
ylabel('R_1 (Hz)')

%% T1 vs depth

for kulma = 1:size(T1_reged.FTDATA,3)
    plotti = new_roi2profile(T1_reged.FTDATA(:,:,kulma),ROI_zones,14);
    depthprofiles(:,kulma) = equisample(plotti,100);
end
[X,Y] = meshgrid(1:100,ang);
figure
hold on
surf(Y,X,depthprofiles')
xlabel('Kulma (aste)')
ylabel('Suhteellinen syvyys ruston pinnalta (%)')
zlabel('T_1 (ms)')

%% Phasemaps


phasemap = zeros(size(ROI_calc.voxels{1}));
h=waitbar(0,['calculating ']);
disced = 0;
for ii = 1:size(T1_reged.FTDATA,1)
    for jj = 1:size(T1_reged.FTDATA,2)
        waitbar(ii/128,h,[num2str(ii), '/' , num2str(128)]);
        
        if (ROI_calc.voxels{1}(ii,jj))
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
                r = squeeze(T1_reged.FTDATA(ii,jj,:)).^-1 - (th(1) + th(2)*(3*cosd(ang + th(3)).^2-1).^2);
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

















