%% Question 1
clear; clc; close all;
delta_t = [0.5, 0.25, 0.125];   % time increment and final time information
global gamma; gamma = 0.75; 
global tau; tau = 1.;

% NR, delta_t = 0.5
tic;
obj1 = NR(@stressR, @stressInc, @HInc, false, delta_t(1));
obj1.solve();
toc;

% NR, delta_t = 0.25
tic;
obj2 = NR(@stressR, @stressInc, @HInc, false, delta_t(2));
obj2.solve();
toc;

% NR, delta_t = 0.125
tic;
obj3 = NR(@stressR, @stressInc, @HInc, false, delta_t(3));
obj3.solve();
toc;

% Modified NR, delta_t = 0.5
tic;
obj4 = NR(@stressR, @stressInc, @HInc, true, delta_t(1));
obj4.solve();
toc;

% Modified NR, delta_t = 0.25
tic;
obj5 = NR(@stressR, @stressInc, @HInc, true, delta_t(2));
obj5.solve();
toc;

% Modified NR, delta_t = 0.125
tic;
obj6 = NR(@stressR, @stressInc, @HInc, true, delta_t(3));
obj6.solve();
toc;

obj_list = [obj1, obj2, obj3, obj4, obj5, obj6];


% Plot the strains
figure
for i = 1:3
   % Cleans up strains so values after NaN don't mess with the plots
   strains = obj_list(i).strain_ev; mask = logical(cumsum(isnan(strains))); strains(mask) = NaN;
   plot(obj_list(i).t_list, strains); hold on    
end
legend('$$\Delta t$$ = 0.5', '$$\Delta t$$ = 0.25', '$$\Delta t$$ = 0.125', 'interpreter', 'latex')
title('Strains of a Non-Linear Viscoelastic Material Computed using Newton-Raphson at Various Time Intervals')
xlabel('Time, $$t$$', 'interpreter', 'latex'); ylabel('Strain, $$\epsilon$$', 'interpreter', 'latex')

figure
for i = 4:6
    % Cleans up strains so values after NaN don't mess with the plots
    strains = obj_list(i).strain_ev; mask = logical(cumsum(isnan(strains))); strains(mask) = NaN;
    plot(obj_list(i).t_list, strains); hold on    
end
legend('$$\Delta t$$ = 0.5', '$$\Delta t$$ = 0.25', '$$\Delta t$$ = 0.125', 'interpreter', 'latex')
title('Strains of a Non-Linear Viscoelastic Material Computed using Modified Newton-Raphson at Various Time Intervals')
xlabel('Time, $$t$$', 'interpreter', 'latex'); ylabel('Strain, $$\epsilon$$', 'interpreter', 'latex')


figure
for i = 1:6
   ratios = log(obj_list(i).residual_ratio);
   if length(ratios) > 10 
        ratios = ratios(1:10);
   end
   plot(1:length(ratios), ratios); hold on
    

end
legend('NR $$\Delta t$$ = 0.5', 'NR $$\Delta t$$ = 0.25', 'NR $$\Delta t$$ = 0.125', ...
       'MNR $$\Delta t$$ = 0.5', 'MNR $$\Delta t$$ = 0.25', 'MNR $$\Delta t$$ = 0.125', 'interpreter', 'latex')
title('Residuals over Iterations Until Convergence at Various Time Intervals')
xlabel('Time, $$t$$', 'interpreter', 'latex'); ylabel('$$\log \left| \frac{R}{R_0} \right|$$', 'interpreter', 'latex')

%%%%%%%%%%% RELEVANT FUNCTIONS %%%%%%%%%%%

% This computes the stress as a function of time
function r = stressR(t)
    if (0 <= t) && (t <= 0.5)
       r = 0.2*t;
    elseif (0.5 < t) && (t <= 1.0)
       r = 0.2*(1-t);
    end
end

% returns stress at next time step
function r = stressInc (cH, dt, ns, cs)
    global gamma; global tau;
    r = ( 1 + gamma * (exp(-dt/(2*tau)) - 1) ) * ( 3 * ns * exp(-10*ns) ) ...
        + 3 * gamma * cH * exp(-dt/tau) ...
        - 3 * gamma * cs * exp(-10*cs) * exp(-dt/(2*tau));
end

% returns history variable at next time step
function r = HInc(Hn, dt, ns, cs)
    global gamma; global tau;
    r = Hn * exp(-dt/tau) + gamma * exp(-dt/(2*tau)) ...
        * ( ns * exp(-10*ns) - cs * exp(-10*cs) );
end