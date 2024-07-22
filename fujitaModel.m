%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjusted from Young et al. 2021

function [shh,svv,shv,svh, modelDepth] = fujitaModel(cof)

  % Constants
  %----------------------------------------------------------------------------------------------------------
  lightSpeed = 299792458;                                                   % Speed of light [m a^-1]
  
  epsilonPerpendicular = 3.15;                                              % Perpendicular component of real part of dielectric permittivity
  deltaEpsilon = 0.034;                                                     % Dielectric anisotropy

  % Radar parameters
  %----------------------------------------------------------------------------------------------------------  
  centerFrequency = 330e6;                                                  % Central frequency [Hz]
  omega = 2*pi*centerFrequency;                                             % Angular frequency

  % model domain
  %----------------------------------------------------------------------------------------------------------  
  azimuthalResolution = 360;                                                % Azimuthal resolution (number of shots)
  numberOfLayers = round(cof.depth(end)-cof.depth(1));                     % Depth resolution (number of numerical layers)
  
  % Define geometry
  %----------------------------------------------------------------------------------------------------------  
  iceThickness = cof.depth(end);                                            % Ice thickness
  startingDepth = cof.depth(1);                                             % starting depth
  depthLimits = [startingDepth,iceThickness];
  
 
  % set up matrices 
  %----------------------------------------------------------------------------------------------------------  

  % Fabric alignment
  %----------------------------------------------------------------------------------------------------------  
  SA(1:4000) = 0;
  phase0g=repmat(SA,size(depthLimits));                                     % Fast axes directions

  
  % Interpolate to numerical domain
  %----------------------------------------------------------------------------------------------------------  
  modelDepth=linspace(startingDepth,iceThickness,numberOfLayers);
  phase0i = phase0g'*pi/180;
  scatteringX = interp1(cof.depth,10.^(cof.rxdBs./10),modelDepth);
  scatteringY = interp1(cof.depth,10.^(cof.rydBs./10),modelDepth);

  
  % Transform fabric eigenvalues to dielectric tensor
  %----------------------------------------------------------------------------------------------------------  
  % Eq. A3 of Fujita et al. (2006)
  % Here eigenvalues are using convention (E3 > E2 > E1)

  modelEpsPerpendicular=epsilonPerpendicular*ones(size(depthLimits));

  Lxi = interp1(cof.depth,cof.exw,modelDepth);
  Lyi = interp1(cof.depth,cof.eyw,modelDepth);

  epsilonX=interp1(depthLimits,modelEpsPerpendicular,modelDepth)+Lxi*deltaEpsilon;
  epsilonY=interp1(depthLimits,modelEpsPerpendicular,modelDepth)+Lyi*deltaEpsilon;
  
  % Calculate average on the layer
  %----------------------------------------------------------------------------------------------------------
  ex=(epsilonX(1:numberOfLayers-1)+epsilonX(2:numberOfLayers))/2;
  ey=(epsilonY(1:numberOfLayers-1)+epsilonY(2:numberOfLayers))/2;
  phase0=(phase0i(1:numberOfLayers-1)+phase0i(2:numberOfLayers))/2;

  Sx=(scatteringX(1:numberOfLayers-1)+scatteringX(2:numberOfLayers))/2;
  Sy=(scatteringY(1:numberOfLayers-1)+scatteringY(2:numberOfLayers))/2;
  
  % Components of T-Matrix
  %----------------------------------------------------------------------------------------------------------
  k0 = 2*pi/(lightSpeed/centerFrequency);                                   % Wave number in a vacuum
  conx = interp1(cof.depth,cof.condx,modelDepth);
  cony = interp1(cof.depth,cof.condy,modelDepth);
  
  condx=(conx(1:numberOfLayers-1)+conx(2:numberOfLayers))/2;                % conductivity in x
  condy=(cony(1:numberOfLayers-1)+cony(2:numberOfLayers))/2;                % conductivity in y

  kx = sqrt((1/lightSpeed^2).*ex.*(omega)^2+1i*4*pi*1e-7.*condx*2*pi*centerFrequency);  % wavenumber x
  ky = sqrt((1/lightSpeed^2).*ey.*(omega)^2+1i*4*pi*1e-7.*condy*2*pi*centerFrequency);  % wavenumber y

  dz=modelDepth(2:end)-modelDepth(1:end-1);
  transmissionXi = exp(-1i.*k0.*dz + 1i.*kx.*dz);                           % Eq. 6a of Fujita et al. (2006)
  transmissionYi = exp(-1i.*k0.*dz + 1i.*ky.*dz);                           % Eq. 6a of Fujita et al. (2006)

  ne = (numberOfLayers-1);                                                  % Number of elements
  R = zeros(2,2,ne);                                                        % Rotation Matrix
  T = zeros(2,2,ne);                                                        % Transmission matrix
  S = zeros(2,2,ne);    

  % Set domain
  %----------------------------------------------------------------------------------------------------------
  x = 0.00001:2*pi/(azimuthalResolution):2*pi;                              % start with 0.00001 to avoid inf.

  for j=1:360

      % Establish matrices
      % -----------------------------------------------------------------------------------------------------
      for i=1:ne
          % Rotation matrix
          % Eq. 10 of Fujita et al. (2006)
          %----------------------------------------------------------------
          R(:,:,i)=[cos(phase0(i)),-sin(phase0(i));sin(phase0(i)),cos(phase0(i))];

          % Scattering matrix
          % Eq. 8 of Fujita et al. (2006)
          %----------------------------------------------------------------
          S(:,:,i)=[Sx(i) 0; 0 Sy(i)];
          % Rotate Scattering
          % Eq. 11 of Fujita et al. (2006)
          S(:,:,i)=R(:,:,i)*S(:,:,i)*R(:,:,i)';

          % Transmission
          % Eq. 5 of Fujita et al. (2006)
          %----------------------------------------------------------------
          T(:,:,i)=[transmissionXi(i) 0; 0 transmissionYi(i)];
          % Rotate transmission
          % Eq. 9 of Fujita et al. (2006)
          T(:,:,i)=R(:,:,i)*T(:,:,i)*R(:,:,i)';
      end

      % Calculate cumulative Transmission (product summation of Eq. 9)
      % cT(n)=Product{R' T R}i=1,n
      %--------------------------------------------------------------------
      cT=T;
      for i=2:ne
          cT(:,:,i)=cT(:,:,i-1)*cT(:,:,i);
      end

      Rt=[cos(x(j)),sin(x(j));-sin(x(j)),cos(x(j))];
      Et=Rt*[1;0];                                                          % transmiter in the frame of reference; [E_PT E_OT] 

      % Eq. 12 of Fujita et al. (2006)
      % -----------------------------------------------------------------------------------------------------

      for i=1:ne
          Er(:,i)=cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Et;                          % Receiver in the frame of reference
          Erant(:,i)=transpose(Rt)*Er(:,i);                                 % Receiver in the frame of antennas
          Erp(i)=Erant(1,i);
          Ero(i)=Erant(2,i);

          % Or reconstruct signal for any azimuthal angle
          %----------------------------------------------------------------
          Err(:,:,i) = cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Rt;
          Errant(:,:,i) = transpose(Rt)*Err(:,:,i);
          Shhn(i) = Errant(1,1,i);
          Svhn(i) = Errant(1,2,i);
          Shvn(i) = Errant(2,1,i);
          Svvn(i) = Errant(2,2,i); 
      end

      ErPar(:,j) = Erp;
      ErPer(:,j) = Ero;

      shh(:,j) = Shhn;
      svh(:,j) = Svhn;
      shv(:,j) = Shvn;
      svv(:,j) = Svvn;
  end
end