clear; clc; close all

global E; E = 79E3;             % in MPa
global K; K = 600;              % in MPa
global H; H = 400;              % in MPa
global sigma_Y; sigma_Y = 200;  % in MPa 
global nu; nu = 0.3;
global mu; mu = 30.38E3;        % in MPa

% Initial variables
delta_t = 0.01; end_time = 4;
times = 0:delta_t:end_time;

eps_p_hist =  zeros(length(times)+1,3,3);
q_hist =  zeros(length(times)+1,3,3);
alpha_hist =  zeros(length(times)+1,1);

stress_hist = zeros(length(times)+1,3,3);
strain_hist = zeros(length(times)+1,3,3);

% Perform return mapping algorithm
for i = 1:length(times)
    [strain_hist(i+1,:,:), stress_hist(i+1,:,:)] = nextStep(times(i));  % temporarily set trial as next stress
    next_dev_tr = squeeze(stress_hist(i+1,:,:)) - trace( squeeze(stress_hist(i+1,:,:)) )/3. .* eye(3);
    next_ksi_tr = next_dev_tr - squeeze(q_hist(i,:,:));
    
    % Check if yielding function is < 0 or not
    next_yF_tr = yieldingFunction( next_ksi_tr, alpha_hist(i) );
    
    if next_yF_tr < 0
        continue;
    else
        next_n = next_ksi_tr/mag(next_ksi_tr);
        delta_gamma = next_yF_tr/(K+H+2*mu);
                
        stress_hist(i+1,:,:) = squeeze(stress_hist(i+1,:,:)) - 2*mu*delta_gamma*next_n;    % Updates stress at next time step  
        eps_p_hist(i+1,:,:) = delta_gamma*next_n + squeeze(eps_p_hist(i,:,:));     % Updates plastic strain at next time step
                
        alpha_hist(i+1) = alpha_hist(i) + delta_gamma;
                
        q_hist(i+1,:,:) = squeeze( q_hist(i,:,:) ) + H*delta_gamma*next_n;
        
    end
    
end

plot(stress_hist(:,2,2))


%%%%%%%%%%% RELEVANT FUNCTIONS %%%%%%%%%%%

% Takes in strain, imposes plane strain conditions, and outputs change in
% sigma
function [r, s] = nextStep(t)
    global E; global nu;
    r = zeros(3); s = zeros(3);
    
    r(2,2) = 0.01*sin(2*pi*t);
    r(1,1) = -nu/(1-nu) * r(2,2);
    
    s(2,2) = E/(1-nu^2) * r(2,2);
    s(3,3) = E*nu/(1-nu^2) * r(2,2);
    
end

function r = yieldingFunction(ksi, alpha)
    global sigma_Y; global K; 
    r = mag(ksi) - ( sqrt(2/3)*sigma_Y + K*alpha );
end

function r = mag(A)
    r = sqrt(sum(A.*A, 'all'));
end




%

