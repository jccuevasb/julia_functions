% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "VfPert" accepts the vortical fissures (VFs) and pertubs them with a "Uniform" or 
% "Gaussian" distribution in the logarithmic domain
% 
% INPUTS:
% ypVF: location of original VFs
% UpVF: velocity of original VFs
% distribution: the distribution of perturbation, Gaussian or Uniform
% percentage: the range of disturbing the VF is    percentage*dlo < y+ < percentage*dyp, 
% where dlo and dup are the gap between each VF and its lower and upper neighbors, respectively
% nIter: number of iterations 
% 
% OUTPUTS:
% ypVFpertsort: sorted location of perturbated VFs
% UPVFpert: velocity of perturbated VFs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ypVfPertSort, UpVfPert] = VfPert(ypVf, UpVf, distribution, percentage, nIter)

nPoints = length(ypVf);
logypVf = log(ypVf);

logypVfpert = zeros(nPoints, nIter);
logypVFpertsort = zeros(nPoints, nIter);
UpVfPert = zeros(nPoints, nIter);

% do not pertub the 1st and the last fissure
logypVfpert(1,:) = logypVf(1);
logypVfpert(end,:) = logypVf(end);

for jj = 1:nIter
    for ii = 2:nPoints - 1
        
        % do not allow perturb VFs go below the lowest and above the highest VF
        while logypVfpert(ii,jj) < logypVf(1) || logypVfpert(ii,jj) > logypVf(end)
            
            % gap between current VF and its neighbours
            logdup = logypVf(ii + 1) - logypVf(ii);
            logdlo = logypVf(ii - 1) - logypVf(ii);
            
            % uniform distribution
            if strcmp(distribution,'Uniform') == 1
                logd = percentage*logdup*(2*(rand(1) - 0.5));
                logypVfpert(ii,jj) = logypVf(ii) + logd;
             
            % Gaussian distribution
            else
                lograndup = percentage*logdup*normrnd(0,1/3);
                lograndlo = percentage*logdlo*normrnd(0,1/3);
                
                % choose whether to move up or down
                temp = rand(1);
                logd = 0;
                while temp == 0;
                    temp = rand(1);
                end
                
                if temp > 0
                    logd = lograndup;
                elseif temp < 0
                    logd = lograndlo;
                end
                logypVfpert(ii,jj) = logypVf(ii)+ logd;
            end
            
        end
    end
    
    % sort veloicties based on location of VFs
    temp = [logypVfpert(:,jj) UpVf];
    tempsort = sortrows(temp,1);
    logypVFpertsort(:,jj) = tempsort(:,1);
    UpVfPert(:,jj) = tempsort(:,2);
end

% transfer back to linear domain 
ypVfPertSort = round(exp(logypVFpertsort));
