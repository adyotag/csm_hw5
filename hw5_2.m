clear; clc; close all

global E; E = 79E3;             % in MPa
global K; K = 600;              % in MPa
global H; H = 400;              % in MPa
global sigma_Y; sigma_Y = 200;  % in MPa 
global nu; nu = 0.3;
global mu; mu = 30.38E3;        % in MPa

% Initial variables
delta_t = 0.001; end_time = 4.;
global times; times = 0:delta_t:end_time;

eps_p_hist =  zeros(length(times)+1,3,3);
q_hist =  zeros(length(times)+1,3,3);
alpha_hist =  zeros(length(times)+1,1);
stress_hist = zeros(length(times)+1,3,3);
strain_hist = zeros(length(times)+1,3,3);

% Perform return mapping algorithm
for i = 1:length(times)
    strain_hist(i+1,:,:) = nextStep(i);  % obtain strain at next step
    delta_strain = squeeze(strain_hist(i+1,:,:) - strain_hist(i,:,:));
    delta_stress = zeros(3); delta_stress(2,2) = E/(1-nu^2) * delta_strain(2,2);
    delta_stress(3,3) = E*nu/(1-nu^2) * delta_strain(2,2);
    stress_hist(i+1,:,:) = squeeze(stress_hist(i,:,:)) + delta_stress; % Obtain trial stress state
        
    next_dev_tr = squeeze(stress_hist(i+1,:,:)) - trace( squeeze(stress_hist(i+1,:,:)) )/3. .* eye(3);
    next_ksi_tr = next_dev_tr - squeeze(q_hist(i,:,:));
    
    next_yF_tr = yieldingFunction( next_ksi_tr, alpha_hist(i) );
        
    % Check if yielding function is < 0 or not
    if next_yF_tr < 0       % Elastic case
        eps_p_hist(i+1,:,:) = eps_p_hist(i,:,:);
        alpha_hist(i+1) = alpha_hist(i);
        
    else        % Plastic case
        next_n = next_ksi_tr/mag(next_ksi_tr);
        delta_gamma = next_yF_tr/(K+H+2*mu);
                
        stress_hist(i+1,:,:) = squeeze(stress_hist(i+1,:,:)) - 2*mu*delta_gamma*next_n;    % Updates stress at next time step  
        eps_p_hist(i+1,:,:) = delta_gamma*next_n + squeeze(eps_p_hist(i,:,:));             % Updates plastic strain at next time step
        alpha_hist(i+1) = alpha_hist(i) + delta_gamma;                                     % Updates alpha at next time step          
        q_hist(i+1,:,:) = squeeze( q_hist(i,:,:) ) + H*delta_gamma*next_n;                 % Updates back stress at next time step
        
    end
    

end

% Plots
figure;
plot(strain_hist(:,1,1),stress_hist(:,1,1)); hold on
plot(strain_hist(:,2,2),stress_hist(:,2,2));
plot(strain_hist(:,3,3),stress_hist(:,3,3));
legend('$$\sigma_{11}$$', '$$\sigma_{22}$$', '$$\sigma_{33}$$', 'interpreter', 'latex');
title('Cyclic Stress-Strain Response under Plane Strain Conditions');
ylabel('Stress, $$\sigma$$ (MPa)', 'interpreter', 'latex');
xlabel('Strain, $$\epsilon$$', 'interpreter', 'latex');

figure;
plot(strain_hist(:,2,2),strain_hist(:,1,1));
title('Lateral Strain against Axial Strain');
ylabel('Lateral Strain, $$\epsilon_{11}$$', 'interpreter', 'latex');
xlabel('Axial Strain, $$\epsilon_{22}$$', 'interpreter', 'latex');

figure;
plot(strain_hist(:,2,2), eps_p_hist(:,2,2));
title('Axial Plastic Strain against Total Axial Strain');
ylabel('Strain, $$\epsilon_{22}^p$$', 'interpreter', 'latex');
xlabel('Total Axial Strain, $$\epsilon_{22}$$', 'interpreter', 'latex');

figure;
plot(strain_hist(:,2,2), alpha_hist);
title('Cumulative Plastic Strain against Total Axial Strain');
ylabel('Cumulative Plastic Strain, $$\alpha$$', 'interpreter', 'latex');
xlabel('Total Axial Strain, $$\epsilon_{22}$$', 'interpreter', 'latex');


%%%%%%%%%%% RELEVANT FUNCTIONS %%%%%%%%%%%

% Imposes plane strain conditions, and outputs strain at next time step
function r = nextStep(i) % Strain, stress
    global nu; global times;
    r = zeros(3); t = times(i);   
    r(2,2) = 0.01*sin(2*pi*t);  % Computes new strain
    r(1,1) = -nu/(1-nu) * r(2,2);
end

% Evaluates the yielding function
function r = yieldingFunction(ksi, alpha)
    global sigma_Y; global K; 
    r = mag(ksi) - ( sqrt(2/3)*sigma_Y + K*alpha );
end

% Returns the magnitude of a tensor
function r = mag(A)
    r = sqrt(sum(A.*A, 'all'));
end


%