clear all
close all

%path2file = '/media/tamara/tamaraPhd/Seagate/EGRIPradarData/';
path2file = '/run/user/1000/gvfs/dav:host=io.erda.dk,ssl=true/ice_radar_data/greenland/greenland_2022/UWB_processed/';
% 
% % 1) plot full profile 20220706
% 
% 
VVlat1 = h5read([path2file,'20220705_103058_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lat');
VVlon1 = h5read([path2file,'20220705_103058_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lon');

VVlat2 = h5read([path2file,'20220705_110610_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lat');
VVlon2 = h5read([path2file,'20220705_110610_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lon');

VVlat3 = h5read([path2file,'20220705_135130_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
VVlon3 = h5read([path2file,'20220705_135130_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');

% VVlat4 = h5read([path2file,'20220705_170928_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
% VVlon4 = h5read([path2file,'20220705_170928_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');
% 
% VVlat5 = h5read([path2file,'20220705_173406_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
% VVlon5 = h5read([path2file,'20220705_173406_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');
% 
% VVlat6 = h5read([path2file,'20220705_174108_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
% VVlon6 = h5read([path2file,'20220705_174108_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');
% 
% VVlat7 = h5read([path2file,'20220705_175514_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
% VVlon7 = h5read([path2file,'20220705_175514_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');
% 
VVlat8 = h5read([path2file,'20220705_180134_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lat');
VVlon8 = h5read([path2file,'20220705_180134_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/lon');

plot(VVlon1(1),VVlat1(1),'*');hold on
plot(VVlon2(1),VVlat2(1),'*');hold on
plot(VVlon3(1),VVlat3(1),'*');hold on
%plot(VVlon4,VVlat4,'*');hold on
%plot(VVlon5,VVlat5,'*');hold on
%plot(VVlon6,VVlat6,'*');hold on
%plot(VVlon7,VVlat7,'*');hold on
plot(VVlon8(1),VVlat8(1),'*');hold on
legend('1','2','3','4');

lat = [VVlat1;VVlat2;VVlat3;VVlat8];
lon = [VVlon1;VVlon2;VVlon3;VVlon8];

wgs84 = wgs84Ellipsoid;
dist = [0;cumsum(distance(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end),wgs84))]; % distance along profile
interv = 0:5000:dist(end);
[a, Idx] = min(abs(dist-interv));
 
% 

%% part1
close all
%clear all
path2file = '/run/user/1000/gvfs/dav:host=io.erda.dk,ssl=true/ice_radar_data/greenland/greenland_2022/UWB_processed/';
profile{1} = '20220705_103058';
profile{2} = '20220705_110610';
profile{3} = '20220705_135130';
profile{4} = '20220705_180134';

v1 = length(VVlat1);
v2 = length(VVlat1)+length(VVlat2);
v3 = length(VVlat1)+length(VVlat2)+length(VVlat3);

Indices{1} = Idx(Idx<=length(VVlat1));
p2 = Idx(Idx<=length(VVlat2)+length(VVlat1));
Indices{2} = p2(p2>max(Indices{1}));
p3 = Idx(Idx<=length(VVlat3)+length(VVlat2)+length(VVlat1));
Indices{3} = p3(p3>max(Indices{2}));
p4 = Idx(Idx>=length(VVlat3)+length(VVlat2)+length(VVlat1));
Indices{4} = p4;


Indices{2} = Indices{2}-v1;
Indices{3} = Indices{3}-v2;
Indices{4} = Indices{4}-v3;

for pcount = 1:length(profile)
%   if pcount>3
%     break
%   end
  clear drivingDir drivingDirection
  profilename = profile{pcount}
  idx = Indices{pcount}

  VVlat = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lat');
  VVlon = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/lon');

  HHtime=  h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/_time');
  VVtime=  h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/_time');

  HH_chirp = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/Chirps');
  VV_chirp = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/Chirps');

  HV_chirp = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_1/Chirps');
  VH_chirp = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_3/Chirps');

  VVdist = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_0/distance');
  HHdist = h5read([path2file,profilename,'_UWB_Greenland_2022_comb.h5'],'/processed/combined/Integrator_2/distance');

  % wgs84 = wgs84Ellipsoid;
  % dist = [0;cumsum(distance(VVlat(1:end-1),VVlon(1:end-1),VVlat(2:end),VVlon(2:end),wgs84))]; % distance along profile
  % 
  % % find indeces at 5 km interval
  % interv = 0:5000:dist(end);
  % [a, idx] = min(abs(dist-interv));
  drivingDir = azimuth(VVlat(1:end-1),VVlon(1:end-1),VVlat(2:end),VVlon(2:end),wgs84);  % azimuth, clockwise from North
  %%
  close all
  clear S


  %genvarname('prof',profilename).idx= idx;
  %%genvarname(['prof',profilename]).lat = VVlat(idx);
  %genvarname(['prof',profilename]).lon = VVlon(idx);
  %prof20220706_163733.dPhh = struct;
  %prof20220706_163733.dPhv = struct;
  eiPar = 3.15; % Parallel component of real part of dielectric permittivity
  eiPer = 3.1840; % Perpendicular component of real part of dielectric permittivity
  c = 299792458; % Speed of light [m a^-1]
  fc = 330e6; % Central frequency [Hz]
  omega = 2*pi*fc; % Angular frequency
  ei = mean([eiPar eiPer]);
  dei = 0.034; %eiPar-eiPer; % Dielectric anisotropy
  X = 360; % Azimuthal resolution (number of shots)

  x = 0:2*pi/(X-1):2*pi;
  xdeg = x * 180/pi;

  theta = 1;
  dzPowerNorm = 50;
  z = HHtime.*1e6.*168./2;
  n=14190;
  dZ = 0.336;
  DptAvgC = dzPowerNorm; 
  stwnd = DptAvgC/dZ;
  S = struct('idx',idx,'lat',VVlat(idx),'lon',VVlon(idx));
  S.z= z(1:n);

  S.time = HHtime(1:n);
  F = figure('units','normalized','outerposition',[0 0 1 1])
  trav= 10; % averaging over 2*trav traces
  for I = 1:length(idx)
  if I==14
    bla=0
  end
    i = idx(I);
    if i>2*trav
      if i>length(HH_chirp.r(:,1))-trav
         
         drivingDirection(I) = median(drivingDir(i-2*trav:i));

         shhs = HH_chirp.r(i-2*trav:i,:)' + 1i*HH_chirp.i(i-2*trav:i,:)';
         svvs = VV_chirp.r(i-2*trav:i,:)' + 1i*VV_chirp.i(i-2*trav:i,:)';  
         svhs = VH_chirp.r(i-2*trav:i,:)' + 1i*VH_chirp.i(i-2*trav:i,:)';
         shvs = HV_chirp.r(i-2*trav:i,:)' + 1i*HV_chirp.i(i-2*trav:i,:)';
      else
      drivingDirection(I) = median(drivingDir(i-trav:i+trav));  
      shhs = HH_chirp.r(i-trav:i+trav,:)' + 1i*HH_chirp.i(i-trav:i+trav,:)';
      svvs = VV_chirp.r(i-trav:i+trav,:)' + 1i*VV_chirp.i(i-trav:i+trav,:)';  
      svhs = VH_chirp.r(i-trav:i+trav,:)' + 1i*VH_chirp.i(i-trav:i+trav,:)';
      shvs = HV_chirp.r(i-trav:i+trav,:)' + 1i*HV_chirp.i(i-trav:i+trav,:)';
      end
%       shhs = HH_chirp.r(i,:)' + 1i*HH_chirp.i(i,:)';
%       svvs = VV_chirp.r(i,:)' + 1i*VV_chirp.i(i,:)';  
%       svhs = VH_chirp.r(i,:)' + 1i*VH_chirp.i(i,:)';
%       shvs = HV_chirp.r(i,:)' + 1i*HV_chirp.i(i,:)';
    else
      drivingDirection(I) = median(drivingDir(i:i+2*trav));
      shhs = HH_chirp.r(1:i+2*trav,:)' + 1i*HH_chirp.i(1:i+2*trav,:)';
      svvs = VV_chirp.r(1:i+2*trav,:)' + 1i*VV_chirp.i(1:i+2*trav,:)';  
      svhs = VH_chirp.r(1:i+2*trav,:)' + 1i*VH_chirp.i(1:i+2*trav,:)';
      shvs = HV_chirp.r(1:i+2*trav,:)' + 1i*HV_chirp.i(1:i+2*trav,:)';
%       shhs = HH_chirp.r(i,:)' + 1i*HH_chirp.i(i,:)';
%       svvs = VV_chirp.r(i,:)' + 1i*VV_chirp.i(i,:)';  
%       svhs = VH_chirp.r(i,:)' + 1i*VH_chirp.i(i,:)';
%       shvs = HV_chirp.r(i,:)' + 1i*HV_chirp.i(i,:)';
    end

    for N = 1:n
      Sn =  [mean(shhs(N,:)) mean(svhs(N,:));mean(shvs(N,:)), mean(svvs(N,:))] ;
      t=1;
      for gamma = 1:360
        thetanew(t) = theta - gamma;
        R = [cosd(thetanew(t)) -sind(thetanew(t)); sind(thetanew(t)), cosd(thetanew(t))];
        Srot = R*Sn*R';
        VV(N,t)= Srot(2,2);
        HH(N,t) = Srot(1,1);
        HV(N,t) = Srot(1,2);
        VH(N,t) = Srot(2,1);
        t=t+1;
      end
    end

    PrParsynt=20*log10(abs(HH));
    PrParsynt=detrend(PrParsynt','constant')';
    PrParsynt=AverageDepth(PrParsynt,z(1:n),dzPowerNorm);

    PrPersynt=20*log10(abs(HV));
    PrPersynt=detrend(PrPersynt','constant')';
    PrPersynt=AverageDepth(PrPersynt,z(1:n),dzPowerNorm);


    chhvv1s = movsum(HH.*conj(VV),[stwnd],1,'omitnan');
    chhvv2s = sqrt(movsum(abs(HH).^2,[stwnd],1,'omitnan'));
    chhvv3s = sqrt(movsum(abs(VV).^2,[stwnd],1,'omitnan'));
    chhvvs = chhvv1s ./ (chhvv2s .* chhvv3s);
    chhvvs = conj(chhvvs);

    [fxs,fys] = gradient(angle(chhvvs));
    phase_ders = 2*c*sqrt(ei)/(4*pi*fc*dei).*fys;


    subplot(2,length(idx),I)
    imagesc(0:360,z(1:n),PrParsynt)
    caxis([-5 5])
    title('dPhh')
    subplot(2,length(idx),length(idx)+I)
    imagesc(0:360,z(1:n),PrPersynt)
    colormap((brewermap([],"RdBu")));
    caxis([-5 5])
    title('dPhv')

    S.dPhv{I} = PrPersynt;
    S.dPhh{I} = PrParsynt;
    S.chhvv{I} = chhvvs;
    S.phasder{I} = phase_ders;
    
  end
  S.drivingDirection = drivingDirection;
  saveas(F,['synth_',profilename,'_',num2str(2*trav),'trace_av.png'])
  figure()
  plot(VVlon,VVlat,'*k');hold on
  plot(VVlon(idx),VVlat(idx),'r*');

  save(['prof',profilename,'_',num2str(2*trav),'trace_av.mat'],'S')

end

%save('prof20220705_attr.mat',')
%% Function for averaging depth between two layers
function y=AverageDepth(x,z,dz)

for i=1:length(z)
    iok=find(z>z(i)-dz/2&z<z(i)+dz/2);
    y(i,:)=mean(x(iok,:),1);
end

end
