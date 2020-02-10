% multiple slab pumping for multipass architecture
% 
% no doping gradient in this version
% 
% multispectral approach
% 
% Daniel Albach                                         2016/11/03

clear all

% at first read the spectra
data_a = dlmread('JK_CaF300Ka.txt');
data_e = dlmread('JK_CaF300Kf.txt');

sigma.l = 900:0.2:1100; % nm!

sigma.a = interp1(data_a(:,1),data_a(:,2),sigma.l);
sigma.e = interp1(data_e(:,1),data_e(:,2),sigma.l);

clear data_a data_e

% only for this particular data multiply everything with 1e-21
% not necessary for Data by JK
% sigma.a = sigma.a.*1e-20;
% sigma.e = sigma.e.*1e-20;

% and wl with 1e-9
% sigma.l = sigma.l.*1e-9;

% mode selection
mode.both = 1; % pump from both sides

% generate geometry

% slab numbers
slab.n = 4;
% doping of each slab (each can have a doping gradient!)
% counting from left to right
slab.dop(1,1) = 1.9;
slab.dop(2,1) = 1.9;
slab.dop(3,1) = 1.9;
slab.dop(4,1) = 1.9;
% thickness of each slab
slab.thick(1) = 0.5; %cm
slab.thick(2) = 0.5; %cm
slab.thick(3) = 0.5; %cm
slab.thick(4) = 0.5; %cm

% slab distance
slab.dist = 0.15; %cm
slab.distabs = 0.0; %absorption coefficient between the slabs [cm^-1]

% spatial and temporal gridding info about the slabs and dist
steps.t = 100;
steps.slab = 100;
steps.dist = 10;

% z-griding points
grid.t = 0:1/(steps.t-1):1;
grid.slab = 0:1/(steps.slab-1):1;
grid.dist = 0:1/(steps.dist-1):1;

% z-griding volume indicators
grid.slabV = ones(steps.slab-1,1);
grid.distV = zeros(steps.dist-1,1);

% constants
const.c = 3e8;
const.N1per = 2.45e20; %for Yb:CaF2; 1.388e20 for Yb:YAG
const.h = 6.626e-34; %Planck
const.tau = 2.4e-3; %fluorescence life time of Yb:CaF2
const.t = 4e-3; %observation time window

% pump informations
pump1.Itot = 8e3; %W/cm²
pump1.T = 4e-3; %s
pump1.cwl = 935e-9; %m
pump1.bw = 1.5e-9; %m

pump2.Itot = 8e3; %W/cm²
pump2.T = 4e-3; %s
pump2.cwl = 935e-9; %m
pump2.bw = 1.5e-9; %m

% the way to calculate all of this:

% generate the geometry
% dist,slab1,dist,slab2,dist,slab3,dist,slab4,dist
z_mesh = grid.dist*slab.dist;
z_mesh = horzcat(z_mesh,grid.slab(2:end)*slab.thick(1)+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.dist(2:end)*slab.dist+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.slab(2:end)*slab.thick(2)+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.dist(2:end)*slab.dist+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.slab(2:end)*slab.thick(3)+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.dist(2:end)*slab.dist+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.slab(2:end)*slab.thick(4)+max(z_mesh));
z_mesh = horzcat(z_mesh,grid.dist(2:end)*slab.dist+max(z_mesh));

% generate the geometry volumes
% dist,slab1,dist,slab2,dist,slab3,dist,slab4,dist
z_meshV = horzcat(grid.distV',grid.slabV'.*1,grid.distV',grid.slabV'.*2,grid.distV',grid.slabV'.*3,grid.distV',grid.slabV'.*4,grid.distV');
z_meshV = int32(z_meshV);


% this turning is very important!
z_mesh = z_mesh';
z_meshV = z_meshV';

z_meshVdist = zeros(size(z_meshV,1),1);

% generate the distances for the volumes
for isteps=1:size(z_meshV,1)
    z_meshVdist(isteps,1) = z_mesh(isteps+1,1)-z_mesh(isteps,1);
end

% do the doping map
dop  = vertcat(zeros(steps.dist-1,1),ones(steps.slab,1)*slab.dop(1),zeros(steps.dist-2,1),ones(steps.slab,1)*slab.dop(2),zeros(steps.dist-2,1),ones(steps.slab,1)*slab.dop(3),zeros(steps.dist-2,1),ones(steps.slab,1)*slab.dop(4),zeros(steps.dist-1,1));
dopV = vertcat(grid.distV,grid.slabV.*slab.dop(1),grid.distV,grid.slabV.*slab.dop(2),grid.distV,grid.slabV.*slab.dop(3),grid.distV,grid.slabV.*slab.dop(4),grid.distV);

% give a diode spectra
% if you don't want to interpolate again, be careful about the spacing and
% I the wavelength!
pump1.dl = 0.2;
pump1.l = pump1.cwl*1e9-20:pump1.dl:pump1.cwl*1e9+20; % nm!
pump1.l = pump1.l.*1e-9;

% I gaussian pump spectrum
pump1.I = exp(-((pump1.l-pump1.cwl).^2)/(2*(pump1.bw)^2));
pump1.I = pump1.I.*pump1.Itot./trapz(pump1.l,pump1.I); % W/m/cm²

% II the wavelength!
pump2.dl = 0.2; % resolution in[nm]
pump2.l = pump2.cwl*1e9-20:pump2.dl:pump2.cwl*1e9+20; % nm!
pump2.l = pump2.l.*1e-9;

% II gaussian pump spectrum
pump2.I = exp(-((pump2.l-pump2.cwl).^2)/(2*(pump2.bw)^2));
pump2.I = pump2.I.*pump2.Itot./trapz(pump2.l,pump2.I); % W/m/cm²

% now make the temporal pump slices
for itime = 1:steps.t
    pump1.Int(itime,:) = pump1.I(:);
    pump2.Int(itime,:) = pump2.I(:);
end

% define the temporal pump structure
% we expect that there are up to two pumps possible
% one from the "left", the other from the "right"
% both might have a different intensity and spectrum
% they might be recycled (incident from the other side again)

% the right/left is done by using a +1 for "from the left" and -1 "from the
% right"

pump1.Intmode = [1,-1];
pump2.Intmode = [-1,1]; %deactivate this field if you only want the pump from one side
% pump1.Intmode = 1;
% pump2.Intmode = -1;

%% clear up all the unecessary stray variables and prepare them for the
% structs to give to the outside functions

mesh.p = z_mesh;
mesh.V = z_meshV;
mesh.Vdist = z_meshVdist;

clear z_mesh z_meshV z_meshVdist

dop_old = dop;
clear dop

dop.p = dop_old;
dop.V = dopV;

clear dop_old dopV

% clean up finished

%% now call the pumping dn/dt

% prepare the initial "zero time" betas
beta.p = zeros(length(mesh.p),1);
beta.V = zeros(length(mesh.V),1);

beta_store = zeros(length(mesh.p),steps.t);

dt = const.t/(steps.t-1);

for itime=1:steps.t
    % generate the average volume beta distribution
    % in the area, where there is no "dist", calculate the beta
    for isteps=1:size(beta.V,1)
        if mesh.V(isteps) ~= 0
            beta.V(isteps) = (beta.p(isteps) + beta.p(isteps+1))/2;
        else
            beta.V(isteps) = 0;
        end
    end

    if isfield(pump1,'Intmode') == 1
        % pump intensity slice
        pump1.It = pump1.Int(itime,:);

        for imode = 1:size(pump1.Intmode,2)

            [dndt_res(:,imode),Int] = dndt(sigma,const,pump1,mesh,dop,beta,slab,pump1.Intmode(imode));

            pump1.It = Int;

        end
    end

    if isfield(pump2,'Intmode') == 1
        % pump intensity slice
        pump2.It = pump2.Int(itime,:);

        for imode = 1:size(pump2.Intmode,2)

            [dndt_res2(:,imode),Int] = dndt(sigma,const,pump2,mesh,dop,beta,slab,pump2.Intmode(imode));

            pump2.It = Int;

        end
    end

    % now do the integration to solve the whole diff equation for one timestep

    dndt_pump = sum(dndt_res,2)+sum(dndt_res2,2);
    for ipos=1:length(mesh.p)
        if dop.p(ipos,1) ~=0.0
            beta.p(ipos,1) = beta.p(ipos,1)+(dndt_pump(ipos,1)/(dop.p(ipos,1)*const.N1per)-beta.p(ipos,1)/const.tau)*dt;
        end
    end
    
    beta_store(:,itime) = beta.p(:,1);
        
end

%% stored energy density
E_den_stored =  trapz(mesh.p,beta.p.*dop.p*const.N1per)*1.95e-19;

%% small signal gain at 1030nm per pass
sigma_abs = interp1(sigma.l,sigma.a,1030);
sigma_ems = interp1(sigma.l,sigma.e,1030);

gainmap = -(sigma_abs - beta.p.*(sigma_abs+sigma_ems)).*dop.p*const.N1per;
SSG = exp(trapz(mesh.p,gainmap));
