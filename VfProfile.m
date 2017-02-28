% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "VfProfile" accepts the vortical fissures (VFs) and generates a master
% velocity profile
% It also groups overlapping VFs to generate a larger VF with the average velocity
%
% INPUTS:
% ypVf: centroid  of  VFs, each column is one iteration
% UpVF: velocity of VFs, each column is one iteration
% h: thickness of the VF
%
% OUTPUTS:
% ypMaster: y-location of master profile
% UpMaster: velocity of master profile, each column is one iteration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ypMaster, UpMaster] = VfProfile(ypVf, UpVf, h)

[~, nIter] = size(UpVf);

nGrid = ypVf(end, 1) + h/2; % grid points in the BL
ypMaster = [1:1:nGrid]';

% fill out the velocity field. the lower and the upper halves of the vf
% must have linear velocity profile Y = a(X-X0) + Y0, or U = m*(y-y0) + U0
UpMaster = zeros(nGrid, nIter);

for jj = 1:nIter
    
    
    %%%%%%%%%%%%%%%%%% combine overlapping profiles %%%%%%%%%%%%%%%%%%%%%%%
    [ypVfLoNew, ypVfCtrNew, ypVfUpNew, UpVfNew] = VfCombine(ypVf(:,jj),  UpVf(:,jj), h);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    newpoints = length(UpVfNew);
    
    %%%%%%%%%%%%%%% Velocity profile within the lowest VF: %%%%%%%%%%%%%%%%
    % points below centroid are 0, point above centroid are linearly distributed
    for ii = 1
        
        % coordinates of the centeroid
        x0 = ypVfCtrNew(ii);
        y0 = UpVfNew(ii);
        
        % upper half slope
        mup =  0.5*(UpVfNew(ii + 1) - UpVfNew(ii))/(ypVfUpNew(ii) - ypVfCtrNew(ii));
        
        % velocity profile of lower half
        UpMaster(1:ypVfLoNew(1), jj) = 0;
                
        % velocity profile of upper half
        UpMaster(ypVfCtrNew(ii):ypVfUpNew(ii), jj) = mup*(ypMaster(ypVfCtrNew(ii):ypVfUpNew(ii)) - x0) + y0;
        
        % velocity profile of the uniform momentum zone (umz) between vfs
        UpMaster(ypVfUpNew(ii):ypVfLoNew(ii + 1), jj) = 0.5*(UpVfNew(ii) + UpVfNew(ii + 1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%% Velocity profile within the middle VFs: %%%%%%%%%%%%%%%
    for ii = 2:newpoints - 1
        
        % coordinates of the centeroid
        x0 = ypVfCtrNew(ii);
        y0 = UpVfNew(ii);
        
        % lower half slope
        mlo = 0.5*(UpVfNew(ii) - UpVfNew(ii - 1))/(ypVfCtrNew(ii) - ypVfLoNew(ii));
        
        % upper half slope
        mup =  0.5*(UpVfNew(ii + 1) - UpVfNew(ii))/(ypVfUpNew(ii) - ypVfCtrNew(ii));
        
        % velocity profile of lower half
        UpMaster(ypVfLoNew(ii):ypVfCtrNew(ii), jj) = mlo*(ypMaster(ypVfLoNew(ii):ypVfCtrNew(ii)) - x0) + y0;
        
        % velocity profile of upper half
        UpMaster(ypVfCtrNew(ii):ypVfUpNew(ii), jj) = mup*(ypMaster(ypVfCtrNew(ii):ypVfUpNew(ii)) - x0) + y0;
        
        % velocity profile of the uniform momentum zone (umz) between vfs
        UpMaster(ypVfUpNew(ii):ypVfLoNew(ii + 1), jj) = 0.5*(UpVfNew(ii) + UpVfNew(ii + 1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%% Velocity profile within the highest VF: %%%%%%%%%%%%%%%
    for ii = newpoints
        
        % coordinates of the centeroid
        x0 = ypVfCtrNew(ii);
        y0 = UpVfNew(ii);
        
        % lower half slope
        mlo = 0.5*(UpVfNew(ii) - UpVfNew(ii - 1))/(ypVfCtrNew(ii) - ypVfLoNew(ii));
        
        % velocity profile of lower half
        UpMaster(ypVfLoNew(ii):ypVfCtrNew(ii), jj) = mlo*(ypMaster(ypVfLoNew(ii):ypVfCtrNew(ii)) - x0) + y0;
        
        % velocity profile of the uniform momentum zone (umz) between vfs
        UpMaster(ypVfCtrNew(ii):ypVfUpNew(ii), jj) = UpVfNew(ii);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
