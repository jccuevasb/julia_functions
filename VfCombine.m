% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "VfCombine" combines overlapping vortical fissures  (VFs) o generate larger VF 
% with average velocity of all overlapped VFs
% 
% INPUTS:
% ypVfSort: sorted location of the centroid of the VFs
% UpVf: centroid velocity of the Vfs 
% h: thickness of the VF
% 
% OUTPUTS:
% ypVfLoNew, ypVfCtrNew, ypVfup: lower, centroid and upper part of combined VFs
% UpVfNew: velocity of combined VFs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ypVfLoNew, ypVfCtrNew, ypVfup, UpVfNew] = VfCombine(ypVfSort,  UpVf, h)

% number of VFs
nPoints = length(ypVfSort);

% grid points
nGrid = ypVfSort(end) + h/2;

% temporary matrices of velocity and location
ytemp = zeros(nGrid, nPoints);
utemp = zeros(nGrid, nPoints);

%%%%%%%%%%%%%%%%%%%%% Find the overlapping VFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1:nPoints
    
    % make location of each VF non-zero
    ytemp(ypVfSort(jj) - h/2:1:ypVfSort(jj) + h/2, jj) = 1;
    
    % make velocity of each VF non-nan
    utemp(:,jj) = nan;
    utemp(ypVfSort(jj) - h/2:1:ypVfSort(jj) + h/2, jj) = UpVf(jj);
end

% assign 1 to each y-location if there is any VF, leave the rest 0
yall = zeros(nGrid, 1);
for ii = 1:nGrid
    if sum(ytemp(ii,:)) > 0
        yall(ii) = 1;
    end
end


% where the difference between two consecutive y-location value is 0 (both have value of 1) they are part of the same VF.
% if the difference is not zero, they are part of two non-overlapping VFs
ind = find(diff(yall));
if mod(length(ind), 2) ~= 0
    ind = [ind;ypVfSort(end) + h/2];
end

% number of new (including overlapped and combined) VFs
n = 0.5*length(ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% index of the lower, centeroid and upper part of each new vf
ypVfLoNew = zeros(n,1);
ypVfCtrNew = zeros(n,1);
ypVfup = zeros(n,1);

% velocity of each new combined VF: is the average of the velocities of all
% overlapped VFs
UpVfNew = zeros(n,1);

for ss = 1:n
    ypVfLoNew(ss) = ind(2*ss - 1) + 1; % the first index
    ypVfup(ss) = ind(2*ss); % the second index
    ypVfCtrNew(ss) = round(0.5*(ypVfLoNew(ss) + ypVfup(ss)));
    
    % average of all vfs' velcoities
    uu = reshape(utemp(ind(2*ss - 1) + 1: ind(2*ss),:),numel(utemp(ind(2*ss - 1) + 1: ind(2*ss),:)),1);
    UpVfNew(ss) = nanmean(uu);
end
