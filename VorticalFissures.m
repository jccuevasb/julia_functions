%% define parameters 

clc
clear all
close all

fontsize = 24;
fontname = 'Times New Roman';
linewidth = 3;
markersize = 10;

load LeeMoserChannel

h = 4; % VF thickness (in plus units)

phic = 0.5*(1+sqrt(5));
deltaplus = 5200; %\delta^+
Uplusinf = 26.5;

% the initial VF
yplus0 = phic*sqrt(deltaplus); 
Uplus0 = 0.5*Uplusinf + 3.5; % note "+ 3.5"

nPoints = 9; % number of VFs in the boundary layer

% location of fissures
ypVf = zeros(nPoints,1); 
ypVf(1) = round(yplus0);

% velocity of fissures
UpVf = zeros(nPoints,1); 
UpVf(1) = Uplus0;

for ii = 2:nPoints
    ypVf(ii) = round(phic*ypVf(ii - 1));
    UpVf(ii) = UpVf(ii - 1) + phic^2*log(phic);
end

% index of the lower, center and upper of each  vf
ypVfLo = ypVf - h/2;
ypVfCtr = ypVf;
ypVfUp = ypVf + h/2;

%%  generate master profile (combine overlapping VFs)
[ypMaster, UpMaster] = VfProfile(ypVf, UpVf, h);


%%%%%%%%%%%%%%%%%%%%%%%% plot the master profile %%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(ypMaster, UpMaster, 'o-', 'markersize', 0.5*markersize); hold on;
plot(ypVf, UpVf, 'ro',  'markersize', markersize);
plot(LM.yp, LM.Up, 'k', 'linewidth', linewidth)
legend('Master profile','VFs','Lee & Moser', 'location','northwest')
xlabel('$y^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
ylabel('$U^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
set(gca,'fontsize',fontsize,'fontname',fontname,'xscale','log');
set(gcf,'Position',[10 40 1260 900])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perturb location of VFs 

% number of iterations
nIter = 3000;

% the range of disturbing the VF is    percentage*dlo < y+ < percentage*dyp, 
% where dlo and dup are the gap between each VF and its lower and upper neighbors, respectively
percentage = 3;

distribution = 'Gaussian';
% distribution = 'Uniform';

[ypVFpertsort, UpVFpert] = VfPert(ypVf, UpVf, distribution, percentage, nIter);

%% Generate velocity Profiles
%TEST 1
[ypPert, UpPert] = VfProfile(ypVFpertsort, UpVFpert, h);

% average all perturbed velocities
UpPertAvg = mean(UpPert,2);

figure;
plot(ypMaster, UpMaster,'g-','linewidth',linewidth); hold on;
plot(ypPert, UpPertAvg,'bo','markersize',markersize/2,'linewidth',linewidth); 
plot(ypVf, UpVf,'ro--','markersize',10,'linewidth',linewidth);
plot(LM.yp, LM.Up, 'k', 'linewidth', linewidth)

xlim([1e1 1e4])
ylim([5 28])
legend('non-pertubed','pertubed average','VFs','Lee & Moser','location','northwest')
xlabel('$y^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
ylabel('$U^+$','fontsize',fontsize,'fontname',fontname,'interpreter','latex');
title(['h = ' num2str(h) ', ' distribution ' dist., ' num2str(nIter) ' iters, ' num2str(percentage*100) ' %'],'fontsize',fontsize,'fontname',fontname,'fontweight','normal');
set(gca,'fontsize',fontsize,'fontname',fontname,'xscale','log');
set(gcf,'Position',[10 40 1260 900])
saveas(gcf,['h' num2str(h,'%0.2d') '_' distribution num2str(percentage*100) '.png'])
