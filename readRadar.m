function [VVlat, VVlon, time, HHchirp,VVchirp,HVchirp,VHchirp,dist] = readRadar(name);
    

    
    VVlat = h5read(name,'/processed/combined/Integrator_2/lat');
    VVlon = h5read(name,'/processed/combined/Integrator_2/lon');
    
    time=  h5read(name,'/processed/combined/Integrator_0/_time');
    %VVtime=  h5read(name,'/processed/combined/Integrator_2/_time');

    HHchirp = h5read(name,'/processed/combined/Integrator_0/Chirps');
    VVchirp = h5read(name,'/processed/combined/Integrator_2/Chirps');

    HVchirp = h5read(name,'/processed/combined/Integrator_1/Chirps');
    VHchirp = h5read(name,'/processed/combined/Integrator_3/Chirps');

    dist = h5read(name,'/processed/combined/Integrator_0/distance');
    %HHdist = h5read(name,'/processed/combined/Integrator_2/distance');
    
end