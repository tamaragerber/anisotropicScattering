clear all
close all
clear all
close all
format long

eiPer = 3.15; % Parallel component of real part of dielectric permittivity
eiPar = 3.1840; % Perpendicular component of real part of dielectric permittivity
co = 299792458; % Speed of light [m a^-1]
fc = 330e6; % Central frequency [Hz]
omega = 2*pi*fc; % Angular frequency
ei = mean([eiPar eiPer]);
dei = 0.034; %eiPar-eiPer; % Dielectric anisotropy
X = 360; % Azimuthal resolution (number of shots)

x = 0:2*pi/(X-1):2*pi;
xdeg = x * 180/pi;


%%  1) observed cirlce radargram

%read UWB data
HHdatetime=  h5read('20220624_concatenated_combined.h5','/concatenated/HH/datetime');
VVdatetime=  h5read('20220624_concatenated_combined.h5','/concatenated/VV/datetime');

HH_chirp = h5read('20220624_concatenated_combined.h5','/concatenated/HH/Chirps');
VV_chirp = h5read('20220624_concatenated_combined.h5','/concatenated/VV/Chirps');

HV_chirp = h5read('20220624_concatenated_combined.h5','/concatenated/HV/Chirps');
VH_chirp = h5read('20220624_concatenated_combined.h5','/concatenated/VH/Chirps');

HHlat = h5read('20220624_concatenated_combined.h5','/concatenated/HH/lat');
HHlon = h5read('20220624_concatenated_combined.h5','/concatenated/HH/lon');

VVlat = h5read('20220624_concatenated_combined.h5','/concatenated/VV/lat');
VVlon = h5read('20220624_concatenated_combined.h5','/concatenated/VV/lon');

VVdist = h5read('20220624_concatenated_combined.h5','/concatenated/VV/distance');
HHdist = h5read('20220624_concatenated_combined.h5','/concatenated/HH/distance');

HHtime = h5read('20220624_concatenated_combined.h5','/concatenated/HH/_time');
VVtime = h5read('20220624_concatenated_combined.h5','/concatenated/VV/_time');

startx = 2722; % start of circle
azstart = 196;azend=674;
startx = startx+azstart-1;
endx = startx+azend-azstart; % end of circle 
z = HHtime.*1e6.*168./2;

wgs84 = wgs84Ellipsoid;
az = azimuth(HHlat(1:end-1),HHlon(1:end-1),HHlat(2:end),HHlon(2:end),wgs84);
azpos = [az(startx-1)-360;az(startx:endx);az(endx+1)+360];
azposlin = 0:360;

% interpolate radargram 
HH_chirplin = struct;HH_chirplin.r = [];HH_chirplin.i = [];
VV_chirplin = struct;VV_chirplin.r = [];VV_chirplin.i = [];


for i = 1:length(HHtime)
  HH_chirplin.r(:,i) = interp1(azpos,HH_chirp.r(startx-1:endx+1,i),azposlin,'linear','extrap');
  HH_chirplin.i(:,i) = interp1(azpos,HH_chirp.i(startx-1:endx+1,i),azposlin,'linear','extrap');
end

for i = 1:length(VVtime)
  VV_chirplin.r(:,i) = interp1(azpos,VV_chirp.r(startx-1:endx+1,i),azposlin,'linear','extrap');
  VV_chirplin.i(:,i) = interp1(azpos,VV_chirp.i(startx-1:endx+1,i),azposlin,'linear','extrap');
end

for i = 1:length(HHtime)
  HV_chirplin.r(:,i) = interp1(azpos,HV_chirp.r(startx-1:endx+1,i),azposlin,'linear','extrap');
  HV_chirplin.i(:,i) = interp1(azpos,HV_chirp.i(startx-1:endx+1,i),azposlin,'linear','extrap');
end

for i = 1:length(VVtime)
  VH_chirplin.r(:,i) = interp1(azpos,VH_chirp.r(startx-1:endx+1,i),azposlin,'linear','extrap');
  VH_chirplin.i(:,i) = interp1(azpos,VH_chirp.i(startx-1:endx+1,i),azposlin,'linear','extrap');
end

HHdat = 20*log10(sqrt(HH_chirp.r(startx-1:endx+1,:).^2 + HH_chirp.i(startx-1:endx+1,:).^2))';
VVdat = 20*log10(sqrt(VV_chirp.r(startx-1:endx+1,:).^2 + VV_chirp.i(startx-1:endx+1,:).^2))';

HHdatlin = 20*log10(sqrt(HH_chirplin.r.^2 + HH_chirplin.i.^2))';
VVdatlin = 20*log10(sqrt(VV_chirplin.r.^2 + VV_chirplin.i.^2))';

clear chhvv shh svv shv svh top bottom1 bottom2
format long
dzPowerNorm = 50; % moving average depth in m. 

shh = HH_chirplin.r' + 1i*HH_chirplin.i';
svv = VV_chirplin.r' + 1i*VV_chirplin.i';
svh = VH_chirplin.r' + 1i*VH_chirplin.i';
shv = HV_chirplin.r' + 1i*HV_chirplin.i';

shhs = imgaussfilt(HH_chirplin.r',10) + 1i*imgaussfilt(HH_chirplin.i',10);
svvs = imgaussfilt(VV_chirplin.r',10) + 1i*imgaussfilt(VV_chirplin.i',10);
svhs = imgaussfilt(VH_chirplin.r',10) + 1i*imgaussfilt(VH_chirplin.i',10);
shvs = imgaussfilt(HV_chirplin.r',10) + 1i*imgaussfilt(HV_chirplin.i',10);


% HH power anomaly
PrPar=20*log10(abs(shh));
PrPar=detrend(PrPar','constant')';
PrPar=AverageDepth(PrPar,z(2:end),dzPowerNorm);

% HV power anomaly

PrPer=20*log10(abs(svh(:,1:end)));
PrPer=detrend(PrPer','constant')';
PrPer=AverageDepth(PrPer,z(2:end),dzPowerNorm);

Prvv=20*log10(abs(svv));
Prvv=detrend(Prvv','constant')';
Prvv=AverageDepth(Prvv,z(2:end),dzPowerNorm);

Prhv=20*log10(abs(shv));
Prhv=detrend(Prhv','constant')';
Prhv=AverageDepth(Prhv,z(2:end),dzPowerNorm);

dZ = 0.336;
DptAvgC = dzPowerNorm; 
stwnd = DptAvgC/dZ;
chhvv1 = movsum(shh.*conj(svv),[stwnd],1,'omitnan');
chhvv2 = sqrt(movsum(abs(shh).^2,[stwnd],1,'omitnan'));
chhvv3 = sqrt(movsum(abs(svv).^2,[stwnd],1,'omitnan'));
chhvv = chhvv1 ./ (chhvv2 .* chhvv3);
chhvv = conj(chhvv);

[fx,fy] = gradient(angle(chhvv));
phase_der = 2*co*sqrt(ei)/(4*pi*fc*dei).*fy;

 w = gausswin(1700, 10);
 w = w/sum(w);

for i = 1:360
 dphidz =  2*co*sqrt(ei)/(4*pi*fc*dei).*gradient(angle(chhvv(:,i)));
 phase_derS(:,i) = filter(w,1,dphidz);
end
II = imgaussfilt(phase_derS,10);
imagesc(II)
% calculate phase difference with a 1D convolutional derivative on both
% real and imaginary components of the complex coherence. 

for i = 1:360
 dphidz_im =  (imag(chhvv(:,i)));
 dphidz_re =  (real(chhvv(:,i)));
 re = filter(w,1,dphidz_re);
 im = filter(w,1,dphidz_im);
 
 phase_derS_conv(:,i) = 2*co*sqrt(ei)/(4*pi*fc*dei).*gradient(angle(re + 1i*im));
end


% [fx,fy] = gradient(angle(C_hhvv));
% phase_derChhvv = 2*c*sqrt(ei)/(4*pi*fc*dei).*fy;



%% calculate synthetic response from hh, hv, vv and vh
plotcount=0;
figure()
for idx=32
theta = idx+90;
plotcount = plotcount+1;
for N = 1:11310
  Sn =  [shhs(N,idx) svhs(N,idx);shvs(N,idx), svvs(N,idx)] ;
  t=1;
  for gamma = 1:360
    thetanew(t) = theta + gamma;
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
PrParsynt=AverageDepth(PrParsynt,z(1:11310),dzPowerNorm);


PrPersynt=20*log10(abs(HV));
PrPersynt=detrend(PrPersynt','constant')';
PrPersynt=AverageDepth(PrPersynt,z(1:11310),dzPowerNorm);


subplot(3,6,plotcount)
imagesc(xdeg,z(1:11310),PrPersynt) 
%colorbar()
caxis([-10 10])
title(num2str(idx))
%axis([0 360 0 3000])
axis([0 360 885 2678])
end
colormap((brewermap([],"RdBu")));
figure()
imagesc(xdeg,z(1:11310),imgaussfilt(PrPar,1)) 
caxis([-10 10])
title('dPhh measured')
axis([0 360 885 2678])
colormap((brewermap([],"RdBu")));


PrPersynt=20*log10(abs(HV));
PrPersynt=detrend(PrPersynt','constant')';
PrPersynt=AverageDepth(PrPersynt,z(1:11310),dzPowerNorm);

dZ = 0.336;
DptAvgC = dzPowerNorm; 
stwnd = DptAvgC/dZ;
chhvv1s = movsum(HH.*conj(VV),[stwnd],1,'omitnan');
chhvv2s = sqrt(movsum(abs(HH).^2,[stwnd],1,'omitnan'));
chhvv3s = sqrt(movsum(abs(VV).^2,[stwnd],1,'omitnan'));
chhvvs = chhvv1s ./ (chhvv2s .* chhvv3s);
%chhvvs = conj(chhvvs);

[fxs,fys] = gradient(angle(chhvvs));
phase_ders = 2*co*sqrt(ei)/(4*pi*fc*dei).*fys;

%% Fujita Model


% read fabric data
fabricdata = readmatrix('../../../externalData/dlilien_fujitaModel/EGRIP_cAxis_form.csv');
scattering = readmatrix('../../../externalData/dlilien_fujitaModel/scatteringSmooth.csv');
d = fabricdata(:,5);

%[ds,s] = sort(d);
[ds,s180] = unique(d);


H = d(end); % Ice thickness
h= d(1);
SA(1:2000) = 0;
SA(2001:3000) = 0;
SA(3001:4000) = 0;
X = 360; % Azimuthal resolution (number of shots)
Z = round(d(end)-d(1)); % Depth resolution (number of numerical layers)
dzPowerNorm = 50; % Moving average depth in metres

% Wave transmission
eiPer = 3.15; % Parallel component of real part of dielectric permittivity
eiPar = 3.1840; % Perpendicular component of real part of dielectric permittivity

% Constants
temp = -20; % Temperature, in Celsius
e_0 = 8.8541878128e-12; % Electric permittivity in vacuum [F m^-1]
mu_0 = 4e-7*pi; % Magnetic permeability in vacuum [H m^-1]
c = 299792458; % Speed of light [m a^-1]

% Radar parameters
fc = 330e6; % Central frequency [Hz]
omega = 2*pi*fc; % Angular frequency

% Define geometry
zg = [h,H];

% Fabric alignment
phase0g=repmat(SA,size(zg)); % Fast axes directions
%phase0g(1:1154) = 90-33;
%phase0g(1155:end) = 157;

% Wave transmission
ei = mean([eiPar eiPer]);
dei = 0.034; %eiPar-eiPer; % Dielectric anisotropy

% Scattering parameters
Sxg = repmat(1,size(zg));
Syg = Sxg ;%.* A;

%%Prepare model input

% Interpolate to numerical domain
zmod=linspace(h,H,Z);
phase0i = phase0g'*pi/180 + 32*pi/180;
[u,q] = unique(scattering(:,1));
a = interp1(u,scattering(q,2),zmod);a(end)=a(end-1);
%test increased scattering in glacial
%a(1400:end) = a(1400:end)/2;
A = 10.^(a./20+log10(1e-12));

Sxi=interp1(zg,Sxg,zmod).*1e-12;
Syi=interp1(zg,Syg,zmod).*A;


% Transform fabric eigenvalues to dielectric tensor
% Eq. A3 of Fujita et al. (2006)
% Here eigenvalues are using convention (E3 > E2 > E1)
exg=eiPer*ones(size(zg));
eyg=eiPer*ones(size(zg));%+dL'*dei;

Lxi = interp1(ds,fabricdata(s180,24),zmod);
Lyi = interp1(ds,fabricdata(s180,25),zmod);
Lyi(1400:end) = Lyi(1400:end) +0.2;

exi=interp1(zg,exg,zmod)+Lxi*dei;
eyi=interp1(zg,eyg,zmod)+Lyi*dei;% .* firnCi;
% test tronger birefringence in lower part
%eyi(1400:end)  = eyi(1400:end)+0.1;

% Calculate average on the layer
ex=(exi(1:Z-1)+exi(2:Z))/2;
ey=(eyi(1:Z-1)+eyi(2:Z))/2;
phase0=(phase0i(1:Z-1)+phase0i(2:Z))/2;
% test slight rotation in glacial
%phase0(1400:end) = pi/4;
phase0deg = rad2deg(phase0)+32;


Sx=(Sxi(1:Z-1)+Sxi(2:Z))/2;
Sy=(Syi(1:Z-1)+Syi(2:Z))/2;

% Components of T-Matrix
k0 = 2*pi/(c/fc); % Wave number in a vacuum
kx=2*pi*fc/co*sqrt(ex); % Propagation speed along x
ky=2*pi*fc/co*sqrt(ey); % Propagation speed along y
dz=zmod(2:end)-zmod(1:end-1);
T_ix = exp(-1i.*k0.*dz + 1i.*kx.*dz); % Eq. 6a of Fujita et al. (2006)
T_iy = exp(-1i.*k0.*dz + 1i.*ky.*dz); % Eq. 6a of Fujita et al. (2006)

ne=(Z-1); % Number of elements
R=zeros(2,2,ne); %Rotation Matrix
T=zeros(2,2,ne); %transmission matrix

% Set domain
x = 0.00001:2*pi/(X-1):2*pi;
xdeg = x * 180/pi;
dx = x(2:end)-x(1:end-1); 

for j=1:359
    
    % Establish matrices
    for i=1:ne
        % Rotation
        % Eq. 10 of Fujita et al. (2006)
        R(:,:,i)=[cos(phase0(i)),-sin(phase0(i));sin(phase0(i)),cos(phase0(i))];
        
        % Scattering
        % Eq. 8 of Fujita et al. (2006)
        S(:,:,i)=[Sx(i) 0; 0 Sy(i)];
        % Rotate Scattering
        % Eq. 11 of Fujita et al. (2006)
        S(:,:,i)=R(:,:,i)*S(:,:,i)*R(:,:,i)';
        
        % Transmission
        % Eq. 5 of Fujita et al. (2006)
        T(:,:,i)=[T_ix(i) 0; 0 T_iy(i)];
        % Rotate transmission
        % Eq. 9 of Fujita et al. (2006)
        T(:,:,i)=R(:,:,i)*T(:,:,i)*R(:,:,i)';
    end
    
    % Calculate cumulative Transmission (product summation of Eq. 9)
    % cT(n)=Product{R' T R}i=1,n
    cT=T;
    for i=2:ne
        cT(:,:,i)=cT(:,:,i-1)*cT(:,:,i);
    end
    
    Rt=[cos(x(j)),sin(x(j));-sin(x(j)),cos(x(j))];
    Et=Rt*[1;0]; %transmiter in the frame of reference; [E_PT E_OT] 
    
    % Eq. 12 of Fujita et al. (2006)
    for i=1:ne
        Er(:,i)=cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Et; %Receiver in the frame of reference
        Erant(:,i)=transpose(Rt)*Er(:,i); %Receiver in the frame of antennas
        Erp(i)=Erant(1,i);
        Ero(i)=Erant(2,i);
        
        % Or reconstruct signal for any azimuthal angle
        Err(:,:,i) = cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Rt;
        Errant(:,:,i) = transpose(Rt)*Err(:,:,i);
        Shhn(i) = Errant(1,1,i);
        Svhn(i) = Errant(1,2,i);
        Shvn(i) = Errant(2,1,i);
        Svvn(i) = Errant(2,2,i); 
    end
    
    ErPar(:,j) = Erp;
    ErPer(:,j) = Ero;
    
    Shh(:,j) = Shhn;
    Svh(:,j) = Svhn;
    Shv(:,j) = Shvn;
    Svv(:,j) = Svvn;

end

% Obtain power anomaly and phase difference

% Detrend power to get anomaly
PrParmod=20*log10(abs(Shh));
PrParmod=detrend(PrParmod','constant')';
PrParmod=AverageDepth(PrParmod,zmod(2:end),dzPowerNorm);
Svh(:,1)= Svh(:,2);
PrPermod=20*log10(abs(Svh(:,1:end)));
PrPermod=detrend(PrPermod','constant')';
PrPermod=AverageDepth(PrPermod,zmod(2:end),dzPowerNorm);

% Calculate and plot polarimetric difference

stwnd = 1;
chhvv1mod = movsum(Shh.*conj(Svv),[stwnd],1,'omitnan'); % dim=1 for depth averages
chhvv2mod = sqrt(movsum(abs(Shh).^2,[stwnd],1,'omitnan'));
chhvv3mod = sqrt(movsum(abs(Svv).^2,[stwnd],1,'omitnan'));
chhvvmod = chhvv1mod ./ (chhvv2mod .* chhvv3mod);
%chhvvmod = conj(chhvvmod);
%% 
[fxm,fym] = gradient(angle(chhvvmod));
phase_derm = 2*co*sqrt(ei)/(4*pi*fc*dei).*fym;

  

figure()
subplot(1,3,1)
imagesc(phase_ders)
title('synthesized')
caxis([-0.5 0.5])
subplot(1,3,2)
imagesc(phase_derm)
title('modeled')
subplot(1,3,3)
imagesc(phase_der)
title('observed')
colormap(brewermap([],"RdBu"));



%% trace minima in psihhvv synth:
%PrPersynt_az = [PrPersynt(:,360-70:end),PrPersynt(:,1:360-70)];
%phase_ders_az = [phase_ders(:,360-70:end),phase_ders(:,1:360-70)];
PrPersynt_az = PrPersynt;
phase_ders_az = phase_ders;
M = PrPersynt_az(:,90:180);
[a,b] = size(M);
for i = 1:a
 [mival(i), mi(i)] = min(M(i,:));
 %[maval(i), ma(i)] = max(M(i,:));

 if i>1
  if (mi(i)-mi(i-1))<-45
    mi(i) = 90+mi(i);
  elseif mi(i)-mi(i-1)>45
    mi(i) = mi(i)-90;
    
  end
%   if (ma(i)-ma(i-1))<-45
%     ma(i) = 90+ma(i);
%   elseif ma(i)-ma(i-1)>45
%     ma(i) = ma(i)-90;
%     
%   end
 end
%  if mival(i)> - 20
%    mival(i) = NaN;
%    mi(i) = NaN;
%  end
 
%  
%  if maval(i) < 20
%    maval(i) = NaN;
%    ma(i) = NaN;
%  end<
%     
end

figure()
imagesc(PrPersynt_az);hold on
plot(mi+90,1:a,'m.')
plot(mi+180,1:a,'g.')
colormap(brewermap([],"RdBu"));
caxis([-10 10])
legend('vx','vy')
vx = mi+90;
vy =mi+180;


for i=1:length(mi)
  V_x(i) = phase_ders_az(i,vx(i));
  V_y(i) = phase_ders_az(i,vy(i));
end

if mean(V_x) > mean(V_y)
  v1 = vy;
  v2 = vx;
else
  v1 = vx;
  v2 = vy;
end
  
figure()
imagesc(phase_ders_az);hold on
plot(v1,1:a,'m.')
plot(v2,1:a,'g.')
legend('v1','v2')
colormap(brewermap([],"RdBu"));
caxis([-0.2 0.2])


M = phase_ders_az(:,1:180);
[a,b] = size(M);
for i = 1:a
 [mival(i), mi(i)] = min(M(i,:));
 if mival(i)> - 20
   mival(i) = NaN;
   mi(i) = NaN;
 end
 [maval(i), ma(i)] = max(M(i,:));
 
 
 if maval(i) < 20
   maval(i) = NaN;
   ma(i) = NaN;
 end
    
end

mean(mival)
mean(maval)
figure()
imagesc(phase_ders_az);hold on
plot(mi+90,1:a,'m.')
plot(ma+90,1:a,'g.')
legend('v1','v2')

colormap(brewermap([],"RdBu"));
caxis([-0.2 0.2])

%% Plotting
% measured
figure()
subplot(3,4,1)
imagesc(xdeg,HHtime(1:11310).*1e6,imgaussfilt(PrPar(1:11310,:),10)) 
%colorbar()
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])
caxis([-10 10])
xticks([0 90 180 270 360])
title('dP_{HH}')
set(gca,'fontsize',14)
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
subplot(3,4,3)
I = imgaussfilt(angle(chhvv),10);
imagesc(xdeg,HHtime(1:11310).*1e6,angle(chhvv(1:11310,:)));
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])

%colorbar()
title('\phi_{HHVV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
subplot(3,4,2)
imagesc(xdeg,HHtime(1:11310).*1e6,imgaussfilt(PrPer(1:11310,:),10));
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])

%colorbar()
title('dP_{HV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
caxis([-10 10])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
colormap((brewermap([],"RdBu")));
subplot(3,4,4)
imagesc(xdeg,HHtime(1:11310).*1e6,phase_der)
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])

title('\psi_{HHVV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
caxis([-0.2 0.2])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])

% synthetic
subplot(3,4,5)
%azimuth correction
%PrParsynt_az = [PrParsynt(:,360-70:end),PrParsynt(:,1:360-70)];
PrParsynt_az = PrParsynt;
imagesc(xdeg,HHtime(1:11310).*1e6,imgaussfilt(PrParsynt_az,10)) 
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])
%colorbar()
caxis([-10 10])
title('dP_{HH}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
%axis([0 360 0 3000])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
subplot(3,4,7)
I = imgaussfilt(angle(chhvvs),10);
%chhvvs_az = [chhvvs(:,360-70:end),chhvvs(:,1:360-70)];
chhvvs_az = chhvvs;
imagesc(xdeg,HHtime(1:11310).*1e6,angle(chhvvs_az));
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])
%colorbar()
title('\phi_{HHVV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
%axis([0 360 0 3000])
subplot(3,4,6)
%PrPersynt_az = [PrPersynt(:,360-70:end),PrPersynt(:,1:360-70)];
PrPersynt_az = PrPersynt;
imagesc(xdeg,HHtime(1:11310).*1e6,imgaussfilt(PrPersynt_az,10));hold on
plot(v1,HHtime(1:11310).*1e6,'k.')
plot(v2,HHtime(1:11310).*1e6,'k.')
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])
%colorbar()
title('dP_{HV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
caxis([-10 10])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])
subplot(3,4,8)
%phase_ders_az = [phase_ders(:,360-70:end),phase_ders(:,1:360-70)];
phase_ders_az=phase_ders;
imagesc(xdeg,HHtime(1:11310).*1e6,phase_ders_az);hold on
plot(v1,HHtime(1:11310).*1e6,'d','markersize',5,'color',[0,0,0],'markerfacecolor',[0.6,0.6,0.8])
plot(v2,HHtime(1:11310).*1e6,'s','markersize',5, 'color',[0,0,0],'markerfacecolor',[0.9,0.7,0.2])
yticks([11.71, 17.63, 23.55, 29.47])
yticklabels([1000, 1500, 2000, 2500])

caxis([-0.2 0.2])
title('\psi_{HHVV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
axis([0 360 HHtime(1).*1e6 HHtime(11310).*1e6])



% modelled
subplot(3,4,9)
PrParmod_az = [PrParmod(:,32:end),PrParmod(:,1:32)];
imagesc(xdeg,zmod,imgaussfilt(PrParmod,10)) 
c=colorbar('location','southoutside');
c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4) ]
c.Label.String = '[dB]';

caxis([-10 10])
xticks([0 90 180 270 360])
title('dP_{HH}')
set(gca,'fontsize',14)
%axis([0 360 0 3000])
axis([0 360 885 2678])
subplot(3,4,11)
I = imgaussfilt(angle(chhvvs),10);
chhvvmod_az = [chhvvmod(:,32:end),chhvvmod(:,1:32)];
imagesc(xdeg,zmod,angle(chhvvmod_az));
c=colorbar('location','southoutside','Ticks',[-pi,0,pi]);
c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4) ]
%c.Position = [ 0.691891891877689 , 0.05,0.212972979506655, 0.013527575198618] ; % or + [left, bottom, width, height] to place it where you want
caxis([-pi pi])
set(gca,'fontsize',14)
%c.Label.String = '[dB]';

title('\phi_{HHVV}')
xticks([0 90 180 270 360])
axis([0 360 885 2678])
%axis([0 360 0 3000])
subplot(3,4,10)
PrPermod_az = [PrPermod(:,32:end),PrPermod(:,1:32)];
imagesc(xdeg,zmod,imgaussfilt(PrPermod_az,10));
c= colorbar('location','southoutside');
c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4) ]
c.Label.String = '[dB]';
xticks([0 90 180 270 360])
title('dP_{HV}')
caxis([-10 10])
set(gca,'fontsize',14)
axis([0 360 885 2678])

subplot(3,4,12)
phase_derm_az = [phase_derm(:,32:end),phase_derm(:,1:32)];
imagesc(xdeg,zmod,phase_derm_az)
caxis([-1 1])
title('\psi_{HHVV}')
set(gca,'fontsize',14)
xticks([0 90 180 270 360])
axis([0 360 885 2678])
c= colorbar('location','southoutside');
c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4) ]


%% plotting 2
% for each depth window, plot az vs power anomaly for model and for circle
% radargram. OBS probably need to shif the azimuth for modelled data to
% match the circle or vise versa. Then fit a 90 deg sine and 180 deg sine
% for both datasets, fitting the amplitude and the azimuth shift. We don't
% really care about the azimuth shift but about the relative amplitudes
% between 90 deg and 180 deg periodicity. 

%mdl = fittype('a*sin(b*x+c)','indep','x');
%fittedmdl2 = fit(x,y,mdl,'start',[rand(),1/(5/8/(2*pi)),rand()]);

% 180 deg periodicity: a*sin(4*pi*xdeg + b)+c

%  X180 = sin(4*pi*xdeg);
% % 90 deg periodicity
%  X90 = sind(8*pi*xdeg);
%  figure()
%  plot(xdeg,X180);hold on
%  plot(xdeg,X90)

f=figure();
%imagesc(diff)
j= 1;
p=1;
%depthint = [500,1000,1500,2000,2500];
st = 0;
zmod = zmod(1:end-1);
for i = 500:500:2500
  if i == 2000
    bla =0;
  end
  try
    % modelled dPhh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %dup = zmod(j);
  %dlo = zmod(i);
  interv = find(zmod<i); interv=find(zmod(interv)>j);
  dup = zmod(interv(1))
  dlo = zmod(interv(end))

  MM = mean(PrParmod(interv,:));
  M = [MM(32:end),MM(1:31)]
  subplot(5,2,2+(p-1)*2)
  plot(xdeg, M,'k.');hold on  
  
    % fitting sin curves
  
  y = M;
  x = xdeg;
  yu = max(y);
  yl = min(y);
  yr = (yu-yl);                               % Range of ‘y’
  yz = y-yu+(yr/2);
  zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
  per = 2*mean(diff(zx));                     % Estimate period
  ym = mean(y);                               % Estimate offset

  fit180 = @(a,x)  a(1).*(sin(4*pi*x + a(2)))+ a(3); % a*sin(4*pi*xdeg + b)
  fcn180 = @(a) sum((fit180(a,x) - y).^2);                              % Least-Squares cost function
  s180 = fminsearch(fcn180, [yr;  pi/2;  ym])                       % Minimise Least-Squares
%  figure(1)
  l1=plot(xdeg,fit180(s180,xdeg), 'r','linewidth',2);hold on

  fit90 = @(b,x)  b(1).*(sin(8*pi*x + b(2)))+ b(3); % a*sin(4*pi*xdeg + b)
  fcn90 = @(b) sum((fit90(b,x) - y).^2);                              % Least-Squares cost function
  s90 = fminsearch(fcn90, [yr;  pi/2;  ym])                       % Minimise Least-Squares
  %xp = linspace(min(x),max(x));
  l2=plot(xdeg,fit90(s90,xdeg),'b','linewidth',2)
   legend([l1,l2],['A_{180} = ',num2str(round(abs(s180(1)),2))],['A_{90} = ',num2str(round(abs(s90(1)),2))],...
   'location','northoutside','orientation','horizontal')
 legend('boxoff')
%     % modelled dPhv
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   subplot(5,4,4+(p-1)*4)
%   M = mean(PrPermod(interv,:));
%   y = M;
%   x = xdeg;
%   yu = max(y);
%   yl = min(y);
%   yr = (yu-yl);                               % Range of ‘y’
%   yz = y-yu+(yr/2);
%   zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
%   per = 2*mean(diff(zx));                     % Estimate period
%   ym = mean(y);                               % Estimate offset
% 
% %   fit180 = @(a,x)  a(1).*(sin(4*pi*x + a(2)))+ a(3); % a*sin(4*pi*xdeg + b)
% %   fcn180 = @(a) sum((fit180(a,x) - y).^2);                              % Least-Squares cost function
% %   s180 = fminsearch(fcn180, [yr;  pi/2;  ym])                       % Minimise Least-Squares
% % %  figure(1)
% %   plot(xdeg,fit180(s180,xdeg), 'r');hold on
% 
%   fit90 = @(b,x)  b(1).*(sin(8*pi*x + b(2)))+ b(3); % a*sin(4*pi*xdeg + b)
%   fcn90 = @(b) sum((fit90(b,x) - y).^2);                              % Least-Squares cost function
%   s90 = fminsearch(fcn90, [yr;  pi/2;  ym])                       % Minimise Least-Squares
%   %xp = linspace(min(x),max(x));
%  plot(xdeg, M,'k.');hold on 
%   plot(xdeg,fit90(s90,xdeg),'r:')

  axis([0 360 -9 9])
  %title([num2str(round(st)),'-',num2str(round(i)), 'm'])
  %set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  set(gca,'Xticklabel',[])
  box on
  grid on
  %xlabel('Aximuth')
  catch
    dup=NaN;
    dlo = NaN;
  end
  
  
  try
  intervr  = find(z<i); interv=find(z(interv)>j);
  dupr = intervr(1);
  dlor = intervr(end);
  
  
  Mr = mean(PrPar(intervr,1:end-2));
  
    % observed dPhh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(5,2,1+(p-1)*2)
  plot(xdeg, Mr,'k.');hold on   

       % fitting sin curves
  
  y = Mr;
  x = xdeg;
  yu = max(y);
  yl = min(y);
  yr = (yu-yl);                               % Range of ‘y’
  yz = y-yu+(yr/2);
  zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
  per = 2*mean(diff(zx));                     % Estimate period
  ym = mean(y);                               % Estimate offset

  fit180 = @(a,x)  a(1).*(sin(4*pi*x + a(2)))+ a(3); % a*sin(4*pi*xdeg + b)
  fcn180 = @(a) sum((fit180(a,x) - y).^2);                              % Least-Squares cost function
  s180 = fminsearch(fcn180, [yr;  pi/2;  ym])                       % Minimise Least-Squares
%  figure(1)
  l1= plot(xdeg,fit180(s180,xdeg), 'r','linewidth',2);hold on

  fit90 = @(b,x)  b(1).*(sin(8*pi*x + b(2)))+ b(3); % a*sin(4*pi*xdeg + b)
  fcn90 = @(b) sum((fit90(b,x) - y).^2);                              % Least-Squares cost function
  s90 = fminsearch(fcn90, [yr;  pi/2;  ym])                       % Minimise Least-Squares
  %xp = linspace(min(x),max(x));
  
  l2= plot(xdeg,fit90(s90,xdeg),'b','linewidth',2)
  legend([l1,l2],['A_{180} = ',num2str(round(abs(s180(1)),2))],['A_{90} = ',num2str(round(abs(s90(1)),2))],...
   'location','northoutside','orientation','horizontal')
 legend('boxoff')
%   % observed dPhv
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Mr = mean(PrPer(intervr,1:end-1));
%   subplot(5,4,2+(p-1)*4)
%   plot(xdeg, Mr,'k.','markerSize',10);hold on   
% 
%        % fitting sin curves
%   
%   y = Mr;
%   x = xdeg;
%   yu = max(y);
%   yl = min(y);
%   yr = double((yu-yl));                               % Range of ‘y’
%   yz = y-yu+(yr/2);
%   zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
%   per = 2*mean(diff(zx));                     % Estimate period
%   ym = double(mean(y));                               % Estimate offset
% % 
% %   fit180 = @(a,x)  a(1).*(sin(4*pi*x + a(2)))+ a(3); % a*sin(4*pi*xdeg + b)
% %   fcn180 = @(a) sum((fit180(a,x) - y).^2);                              % Least-Squares cost function
% %   s180 = fminsearch(fcn180, [yr;  pi/2;  ym])                       % Minimise Least-Squares
% % %  figure(1)
% %   plot(xdeg,fit180(s180,xdeg), 'b','linewidth',2);hold on
% 
%   fit90 = @(b,x)  b(1).*(sin(8*pi*x + b(2)))+ b(3); % a*sin(4*pi*xdeg + b)
%   fcn90 = @(b) sum((fit90(b,x) - y).^2);                              % Least-Squares cost function
%   s90 = fminsearch(fcn90, [yr;  pi/2;  ym])                       % Minimise Least-Squares
%   %xp = linspace(min(x),max(x));
%   plot(xdeg,fit90(s90,xdeg),'b:','linewidth',2)
  catch
    dupr=NaN;
    dlor =NaN;
    box off
  end
  
%   if isnan(dup)==0
%     title([num2str(round(st)),'-',num2str(round(i)), 'm'])
%   else
%     title([num2str(round(st)),'-',num2str(round(i)), 'm'])
%   end
 % set(gca,'fontsize',14)
  axis([0 360 -9 9])
  xticks([0 90 180 270 360])
  set(gca,'Xticklabel',[])
  
  box on
  grid on
  
%   if p<5
%     set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
%   else
%     xlabel('Aximuth')
%   end
  j= i;
  p=p+1;
  st = i;
  
end
f.Position = [0 0 550 1200];
%%
close all
save('headingFig_input.mat')

%% Function for averaging depth between two layers
function y=AverageDepth(x,z,dz)

for i=1:length(z)
    iok=find(z>z(i)-dz/2&z<z(i)+dz/2);
    y(i,:)=mean(x(iok,:),1);
end

end

