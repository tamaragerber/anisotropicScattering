
% calculate power ratio for anisotropic scattering
% -------------------------------------------------------------------------
function [cof] = calculatePowerRatio(cof)

   meandep = zeros(1,length(cof.ex)-1);
   rw  = zeros(1,length(cof.ex)-1);

    for i = 1:length(cof.ex)-1
        meandep(i) = mean([cof.depth(i), cof.depth(i+1)]);
        rw(i) = abs(((cof.eyw(i) - cof.eyw(i+1)) / (cof.exw(i) - cof.exw(i+1))));
    end
    
    rdB = 10 * log10(rw);                                                   % definition by ershadi 2022
    
    [cof.Zx, cof.Zy, cof.rxdB, cof.rydB] = calculateImpedanceRatio(cof);    % definition by fujita 2000
    
    % clean vectors
    %-----------------------------------------------------------------------
    infIndex = isinf(abs(rdB));
    rdB(infIndex) = [];
    meandep(infIndex) = [];
    cof.rxdB(infIndex) = [];
    cof.rydB(infIndex) = []; 
    
    nanIndex = isnan(abs(rdB));
    rdB(nanIndex) = [];
    meandep(nanIndex) = [];
    cof.rxdB(nanIndex) = [];
    cof.rydB(nanIndex) = []; 
    
    if isempty(cof.rxdB)
      cof.rdB = cof.depth.*0;
      cof.rxdB = cof.depth.*0;
      cof.rydB = cof.depth.*0;
    else
      cof.rdB = interp1(meandep,rdB,cof.depth);
      cof.rxdB = interp1(meandep,cof.rxdB,cof.depth);
      cof.rydB = interp1(meandep,cof.rydB,cof.depth);
    end
    
    cof.rdBs = smooth(cof.rdB,50);
    cof.rxdBs = smooth(cof.rxdB,50);
    cof.rydBs = smooth(cof.rydB,50);
    cof.meandep = meandep;
end

%% additional function

% calculate impedance
% -------------------------------------------------------------------------
function [Zx,Zy, rxdB, rydB] = calculateImpedanceRatio(cof)
    vacuumPermittivity = 8.8541878128e-12;                                  % Electric permittivity in vacuum [F m^-1]
    vacuumPermeability = 4e-7*pi;                                           % Magnetic permeability in vacuum [H m^-1]
    centerFrequency = 330e6;                                                % Central frequency [Hz]
    omega = 2*pi*centerFrequency;                                           % Angular frequency

    eiPer = 3.15;                                                           % Parallel component of real part of dielectric permittivity
    dei = 0.034;          
    
    epsx = eiPer+cof.exw*dei;
    epsy = eiPer+cof.eyw*dei;
    gammax = sqrt(vacuumPermittivity*vacuumPermeability*epsx.*omega^2 + 1i*vacuumPermeability*cof.condx*omega);
    gammay = sqrt(vacuumPermittivity*vacuumPermeability*epsy.*omega^2 + 1i*vacuumPermeability*cof.condy*omega);
    Zx = (1i*vacuumPermeability*omega)./gammax;
    Zy = (1i*vacuumPermeability*omega)./gammay;
    
    rx = zeros(1,length(Zx)-1);
    ry = zeros(1,length(Zx)-1);
    for i =  1:length(Zx)-1
      rx(i) = (Zx(i)-Zx(i+1))/(Zx(i)+Zx(i+1));
      ry(i) = (Zy(i)-Zy(i+1))/(Zy(i)+Zy(i+1));
    end
    
    rxdB = 10*log10(abs(rx));
    rydB = 10*log10(abs(ry));

end
