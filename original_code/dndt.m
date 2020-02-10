% function dndt for the calulation of the dn/dt part for one timestep and
% only one intensity "run"
% 
% this code is developed for readability, not speed optimization
% 
% 
% Daniel Albach                                                 2013/04/12

function [dndt,Int] = dndt(sigma,const,pump,mesh,dopdist,betadist,slab,mode)

z_mesh = mesh.p;
z_meshV = mesh.V;
z_meshVdist = mesh.Vdist;

beta = betadist.p;
beta_vol = betadist.V;

dop = dopdist.p;
dopV = dopdist.V;

% if we want to pump from the other side

if mode ~=1
    beta = flipud(beta);
    beta_vol = flipud(beta_vol);
    dopV = flipud(dopV);
    dop = flipud(dop);
    z_meshVdist = flipud(z_meshVdist);
    z_meshV = flipud(z_meshV);
    z_mesh = flipud(z_mesh);
end


I1 = zeros(length(z_mesh),size(pump.l,2));
dn1 = zeros(length(z_mesh),size(pump.l,2));

for iwl=1:size(pump.l,2)
    
    sigma_abs = interp1(sigma.l,sigma.a,pump.l(iwl)*1e9);
    sigma_ems = interp1(sigma.l,sigma.e,pump.l(iwl)*1e9);
    
    I1(1,iwl) = pump.It(iwl)*pump.dl*1e-9;
    
    dn1(1,iwl) = I1(1,iwl)* ((sigma_abs - beta(1,1)*(sigma_abs+sigma_ems))*dop(1,1)*const.N1per)/(const.h*const.c)*pump.l(iwl);

    for isteps = 2:length(z_mesh)
        if z_meshV(isteps-1) ~= 0
            I1(isteps,iwl) = I1(isteps-1,iwl)*exp(-(sigma_abs - beta_vol(isteps-1,1)*(sigma_abs+sigma_ems))*dopV(isteps-1,1)*z_meshVdist(isteps-1,1)*const.N1per);
        else
            I1(isteps,iwl) = I1(isteps-1,iwl)*exp(-slab.distabs*z_meshVdist(isteps-1,1));
        end
        
        dn1(isteps,iwl) = I1(isteps,iwl)*((sigma_abs - beta(isteps,1)*(sigma_abs+sigma_ems))*dop(isteps,1)*const.N1per)/(const.h*const.c)*pump.l(iwl);
        
    end
end

% calculate the dn/dt mapping and sum it up
% this is the point, where the spectral positioning and resolution is
% kicked out

I1 = I1.*1e9;

if mode ~=1
    dndt = flipud(sum(dn1,2));
    Int = flipud(I1(length(z_mesh),:));
else
    dndt = sum(dn1,2);
    Int = I1(length(z_mesh),:);
end

end