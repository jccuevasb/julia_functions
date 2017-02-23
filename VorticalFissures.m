%% master profile 
 
clear all
% close all

fontsize = 24;
fontname = 'Times New Roman';
linewidth = 3;
markersize = 10;

phic = 0.5*(1+sqrt(5));
deltaplus = 5200; %\delta^+
Uplusinf = 26.5;

yplus0 = phic*sqrt(deltaplus);
Uplus0 = 0.5*Uplusinf + 3.5; % note "+ 3.5"

npoints = 9; % number of fissures in the boundary layer
yplusFiss = zeros(npoints,1); % location of fissures
yplusFiss(1) = round(yplus0);

h = 2; % fissure thickness
%%
UplusFiss = zeros(npoints,1); % velocity of fissures
UplusFiss(1) = Uplus0;

for ii = 2:npoints
    yplusFiss(ii) = round(phic*yplusFiss(ii - 1));
    UplusFiss(ii) = UplusFiss(ii - 1) + phic^2*log(phic);
end
logyplusFiss = log(yplusFiss);


UpFissMaster = zeros(npoints, h+1);
ypFissMaster = zeros(npoints, h+1);
for ii = 1:npoints
    ypFissMaster(ii,:) = [yplusFiss(ii) - h/2:yplusFiss(ii) + h/2];
    UpFissMaster(ii,:) = UplusFiss(ii);
end


% discretize yplus with increments of 1;
ngrid = yplusFiss(end) - yplusFiss(1) + 1;
yplus_nopert = [yplusFiss(1):1:yplusFiss(end)]';
logyplus_nopert = log(yplus_nopert);
Uplus_nopert = zeros(ngrid,1);
Uplus_nopert(1) = UplusFiss(1);


% assign velocities to additional positions
% index of the lower, center and upper of each  vf
yfisslo = yplusFiss - h/2;
yfissctr = yplusFiss;
yfissup = yplusFiss + h/2;

% velocity of each  vf
ufiss = UplusFiss;

% add (0,0) to location and velocity of vfs
newpoints = npoints+1;
yfisslo = [0;yfisslo];
yfissctr = [0;yfissctr];
yfissup = [0;yfissup];

ufiss = [0;ufiss];

% fill out the velocity field. the lower and the upper halves of the vf
% must have linear velocity profile Y = a(X-X0) + Y0, or U = m*(y-y0) + U0

yplus = [1:1:yplus_nopert(end) + h/2]';
Uplus = zeros(yplus_nopert(end) + h/2,1);

for ii = 1
    % velocity profile of the uniform momentum zone (umz) between vfs
    Uplus(yfissup(ii) + 1:yfisslo(ii + 1)) = 0.5*(ufiss(ii) + ufiss(ii + 1));
end

for ii = 2:newpoints - 1
    
    % coordinates of the centeroid
    x0 = yfissctr(ii);
    y0 = ufiss(ii);
    
    % lower half slope
    mlo = 0.5*(ufiss(ii) - ufiss(ii - 1))/(yfissctr(ii) - yfisslo(ii));
    
    % upper half slope
    mup =  0.5*(ufiss(ii + 1) - ufiss(ii))/(yfissup(ii) - yfissctr(ii));
    
    % velocity profile of lower half
    Uplus(yfisslo(ii):yfissctr(ii)) = mlo*(yplus(yfisslo(ii):yfissctr(ii)) - x0) + y0;
    
    % velocity profile of upper half
    Uplus(yfissctr(ii):yfissup(ii)) = mup*(yplus(yfissctr(ii):yfissup(ii)) - x0) + y0;
    
    % velocity profile of the uniform momentum zone (umz) between vfs
    Uplus(yfissup(ii):yfisslo(ii + 1)) = 0.5*(ufiss(ii) + ufiss(ii + 1));
end

for ii = newpoints
    
    % coordinates of the centeroid
    x0 = yfissctr(ii);
    y0 = ufiss(ii);
    
    % lower half slope
    mlo = 0.5*(ufiss(ii) - ufiss(ii - 1))/(yfissctr(ii) - yfisslo(ii));
    
    % velocity profile of lower half
    Uplus(yfisslo(ii):yfissctr(ii)) = mlo*(yplus(yfisslo(ii):yfissctr(ii)) - x0) + y0;
    
    % velocity profile of the uniform momentum zone (umz) between vfs
    Uplus(yfissctr(ii):yfissup(ii)) = ufiss(ii);
end

% figure;
% semilogy(Uplus, yplus,'*-'); hold on;
% semilogy(UplusFiss,yplusFiss,'ro');

%% perturbation (perturb iteratively) 
% number of iterations
nsteps = 3000;
% close all
% the range of disturbing the fissures, is perc*dup and perc*dlo
perc = 3;

% location of perturbed fissures
% do not pertub the 1st and the last fissure
logyplusFisspert = zeros(npoints, nsteps);
logyplusFisspertsort = zeros(npoints, nsteps);
UplusFisspert = zeros(npoints, nsteps);

logyplusFisspert(1,:) = logyplusFiss(1);
logyplusFisspert(end,:) = logyplusFiss(end);

distribution = 'Gaussian';
% distribution = 'Uniform';

jj = 1;
for ii = 2:npoints - 1
    while logyplusFisspert(ii,jj) < logyplusFiss(1) || logyplusFisspert(ii,jj) > logyplusFiss(end)
        logdup = logyplusFiss(ii + 1) - logyplusFiss(ii);
        logdlo = logyplusFiss(ii - 1) - logyplusFiss(ii);
        
        % uniform distribution
        if strcmp(distribution,'Uniform') == 1
            logd = perc*logdup*(2*(rand(1) - 0.5));
            logyplusFisspert(ii,jj) = logyplusFiss(ii) + logd;
            
        else
            % Gaussian distribution
            lograndup = perc*logdup*normrnd(0,1/3);
            lograndlo = perc*logdlo*normrnd(0,1/3);
            
            logyplusFisspert(ii,jj) = logyplusFiss(ii)+ lograndup;
            %                 % choose whether to move up or down
            %                 temp = rand(1);
            %                 logd = 0;
            %                 while temp == 0;
            %                     temp = rand(1);
            %                 end
            %
            %                 if temp > 0
            %                     logd = lograndup;
            %                 elseif temp < 0
            %                     logd = lograndlo;
            %                 end
            %                 logyplusFisspert(ii,jj) = logyplusFiss(ii)+ logd;
        end
        
    end
end
temp = [logyplusFisspert(:,jj) UplusFiss];
tempsort = sortrows(temp,1);
logyplusFisspertsort(:,jj) = tempsort(:,1);
UplusFisspert(:,jj) = tempsort(:,2);

for jj = 2:nsteps
    for ii = 2:npoints - 1
        while logyplusFisspert(ii,jj) < logyplusFiss(1) || logyplusFisspert(ii,jj) > logyplusFiss(end)
            logdup = logyplusFisspertsort(ii + 1, jj - 1) - logyplusFisspertsort(ii, jj - 1);
            logdlo = logyplusFisspertsort(ii - 1, jj - 1) - logyplusFisspertsort(ii, jj - 1);
            
            
            % uniform distribution
            if strcmp(distribution,'Uniform') == 1
                lograndup = perc*logdup*rand(1);
                lograndlo = perc*logdlo*rand(1);
            
            % Gaussian distribution
            else 
                lograndup = perc*logdup*normrnd(0,1/3);
                lograndlo = perc*logdlo*normrnd(0,1/3);
            end
            
            % choose whether to move up or down
            temp = rand(1);
            logd = 0;
            while temp == 0.5;
                temp = rand(1);
            end
            
            if temp > 0.5
                logd = lograndup;
            elseif temp < 0.5
                logd = lograndlo;
            end
            
%             logyplusFisspert(ii,jj) = logyplusFisspertsort(ii,jj - 1) + logd; % iterative pertirbation
            logyplusFisspert(ii,jj) = logyplusFiss(ii) + logd; % perturb the base only
        end          
    end
    temp = [logyplusFisspert(:,jj) UplusFiss];
    tempsort = sortrows(temp,1);
    logyplusFisspertsort(:,jj) = tempsort(:,1);
    UplusFisspert(:,jj) = tempsort(:,2);
end

yplusFisspert = round(exp(logyplusFisspert));
yplusFisspertsort = round(exp(logyplusFisspertsort));

%%
clc

yplus_pert = [1:1:yplus_nopert(end) + h/2]';
Uplus_pert = zeros(yplus_nopert(end) + h/2,nsteps);

% assign velocities to additional positions 
for jj = 1:nsteps
    
    % group overlapping vortical fissures  (vf)
    % put each vf in one column and keep track of location indices
    ytemp = zeros(yplus_pert(end) + h/2, npoints);
    utemp = zeros(yplus_pert(end) + h/2, npoints);
    
    for ii = 1:npoints
        
        % make location each vf non-zero
        ytemp(yplusFisspertsort(ii,jj) - h/2:1:yplusFisspertsort(ii,jj) + h/2, ii) = 1; 
        
        % make velocity of each vf non-nan
        utemp(:,ii) = nan;
        utemp(yplusFisspertsort(ii,jj) - h/2:1:yplusFisspertsort(ii,jj) + h/2, ii) = UplusFisspert(ii,jj);
    end
    
    % average at all y-locations; each location with more than 0 vf will be non-zero
    ym = mean(ytemp,2);
    
    % find non-zero wall locations and make them ones
    z1 = find(ym);
    ym(z1) = 1;
    
    % where the difference between two consecutive y-location value is 0, they are part of the same vf. 
    % if the difference is not zero, they are part of two non-overlapping vfs
    z = find(diff(ym));
    
    % add index of the last vf
    %z = [z;yplus_pert(end) + h/2];
    
    % number of new (including overlapped and combined) vfs
    n = 0.5*numel(z);
    
    % index of the lower, center and upper of each new vf
    yfisslo = zeros(n,1);
    yfissctr = zeros(n,1);
    yfissup = zeros(n,1);
    
    % velocity of each new vf: is the average of the velocities of all
    % overlapped vfs
    ufiss = zeros(n,1);
    
    for ss = 1:n
        yfisslo(ss) = z(2*ss - 1) + 1; % the first index
        yfissup(ss) = z(2*ss); % the second index 
        yfissctr(ss) = floor(0.5*(yfisslo(ss) + yfissup(ss)));
        
        % average of all vfs' velcoities
        uu = reshape(utemp(z(2*ss - 1) + 1: z(2*ss),:),numel(utemp(z(2*ss - 1) + 1: z(2*ss),:)),1);
        ufiss(ss) = nanmean(uu);
    end
    
    % add (0,0) to location and velocity of vfs
    newpoints = n+1;
    yfisslo = [0;yfisslo];
    yfissctr = [0;yfissctr];
    yfissup = [0;yfissup];
    
    ufiss = [0;ufiss];

    % fill out the velocity field. the lower and the upper halves of the vf
    % must have linear velocity profile Y = a(X-X0) + Y0, or U = m*(y-y0) + U0
    
    for ii = 1
        % velocity profile of the uniform momentum zone (umz) between vfs
        Uplus_pert(yfissup(ii) + 1:yfisslo(ii + 1), jj) = 0.5*(ufiss(ii) + ufiss(ii + 1));
    end
    
    for ii = 2:newpoints - 1
        
        % coordinates of the centeroid
        x0 = yfissctr(ii);
        y0 = ufiss(ii);
        
        % lower half slope
        mlo = 0.5*(ufiss(ii) - ufiss(ii - 1))/(yfissctr(ii) - yfisslo(ii));
        
        % upper half slope
        mup =  0.5*(ufiss(ii + 1) - ufiss(ii))/(yfissup(ii) - yfissctr(ii));
        
        % velocity profile of lower half
        Uplus_pert(yfisslo(ii):yfissctr(ii), jj) = mlo*(yplus_pert(yfisslo(ii):yfissctr(ii)) - x0) + y0;
        
        % velocity profile of upper half
        Uplus_pert(yfissctr(ii):yfissup(ii), jj) = mup*(yplus_pert(yfissctr(ii):yfissup(ii)) - x0) + y0;
        
        % velocity profile of the uniform momentum zone (umz) between vfs
        Uplus_pert(yfissup(ii):yfisslo(ii + 1), jj) = 0.5*(ufiss(ii) + ufiss(ii + 1));
    end
    
    for ii = newpoints
        
        % coordinates of the centeroid
        x0 = yfissctr(ii);
        y0 = ufiss(ii);
        
        % lower half slope
        mlo = 0.5*(ufiss(ii) - ufiss(ii - 1))/(yfissctr(ii) - yfisslo(ii));

        % velocity profile of lower half
        Uplus_pert(yfisslo(ii):yfissctr(ii), jj) = mlo*(yplus_pert(yfisslo(ii):yfissctr(ii)) - x0) + y0;
        
        % velocity profile of the uniform momentum zone (umz) between vfs
        Uplus_pert(yfissctr(ii):yfissup(ii), jj) = ufiss(ii);
    end
end

% average all perturbed velocities
Uplus_pert_m = mean(Uplus_pert,2);

%%


% TEST TEST
load yp5200pp
figure;
plot(yplus, Uplus,'g-','linewidth',linewidth); hold on;
plot(yplus_pert, Uplus_pert_m,'bo','markersize',markersize/2,'linewidth',linewidth); 
plot(yplusFiss, UplusFiss,'ro--','markersize',10,'linewidth',linewidth);
plot(yp5200, U5200,'k-','markersize',10,'linewidth',linewidth);

xlim([1e1 1e4])
ylim([5 28])
legend('non-pertubed','pertubed average','VFs','Lee & Moser','location','northwest')
xlabel('$y^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
ylabel('$U^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
title(['h = ' num2str(h) ', ' distribution ' dist., ' num2str(nsteps) ' iters, ' num2str(perc*100) ' %'],'fontsize',fontsize,'fontname',fontname,'fontweight','normal');
set(gca,'fontsize',fontsize,'fontname',fontname,'xscale','log');
set(gcf,'Position',[10 40 1260 900])
saveas(gcf,['h' num2str(h,'%0.2d') '_' distribution num2str(perc*100) '.png'])
