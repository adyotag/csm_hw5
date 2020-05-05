classdef NR < handle
   properties
        init_guess = 0.01;                                            % default inital guess for iterative method
        stressR = NaN; residual_handle = NaN;                       % stressR is the loading function specified in the problem statement, strainR is either the first or second constitutive law
        stressInc = NaN; HInc = NaN; init_strain = 0; init_H = 0;   % stressInc and HInc are functions to solve for stress and history variable at next time step
        modified = false; first_exec = true;                        % if modified is true, then we perform modified Newton Raphson method; first_exec detects whether we performed a first iteration
        t_start = 0; t_end = 1;                                     % start and end times to solve for. For each time step, we iteratively find the solution
        delta_x = 1E-3; delta_t = 0.5;                              % delta_x helps obtain derivatives numerically, and delta_t is what we specify
        strain_ev = NaN; stress_ev = NaN; t_list = NaN;             % strain_ev is the strain we solve for. stress_ev is from the loading stress we specify. Both are functions of time
        first_residual = NaN; residual_ratio = [];                  % saves the first residual; this is needed for looking at convergence in the assignment; residual_ratio stores the ratio ||Rn/R0||
        num_iters = 0;                                              % Tracks the number of iteration that the NR or modified NR went through
        abs_thresh = 1.E-2; rel_thresh = 2.5E-1;           % Absolute and relative thresholds
   end
   methods
       % Constructor
       function self = NR(stressR, stressInc, HInc, modified, delta_t, abs_thresh, rel_thresh, ...
                            init_strain, init_H, t_start, t_end, delta_x,  init_guess)
            
            self.stressR = stressR;                                 % Specify loading function
            self.stressInc = stressInc; self.HInc = HInc;           % Specify nonlinear viscoelasticity incrementation functions
            self.modified = modified;                               % Choose classic or modified NR
            self.delta_t = delta_t;                                 % delta t of choice
            
            if nargin == 7
                self.abs_thresh = abs_thresh;                       % Set init_guess, absolute and relative thresholds if specified
                self.rel_thresh = rel_thresh;
            end
            
            if nargin > 7
                self.t_start = t_start;                                 % Specify time information and delta_x for numerical derivative
                self.t_end = t_end;
                self.delta_x = delta_x;
                self.init_guess = init_guess;                           % Specify initial guess
                self.init_strain = init_strain; self.init_H = init_H;    % Specifies initial strain and history variables
            end
            
            self.t_list = self.t_start:self.delta_t:self.t_end;             % list of all times we need to solve for strains
            self.stress_ev = arrayfun(@self.stressR, self.t_list);          % vectorized computation of stress at each time
            self.strain_ev = zeros(size(self.stress_ev));                   % allocate memory for strain 
       end
       
       % Returns stress response
       function r = getStress(self)
          r = self.stress_ev;
       end
       
       % Returns strain response
       function r = getStrain(self)
           r = self.strain_ev;
       end
       
       % Returns residual ratios
       function r = getResiduals(self)
          r = self.residual_ratio; 
       end

       % Computes strain at each time. For each time, we use the (modified)
       % Newton Raphson method to numerically solve for strain at the next
       % time step. Using this information, we compute the history variable
       % at the next time step. We then set out initial guess to the strain
       % at the next time step, and continue to solve for strain.
       function solve(self)
           cH = self.init_H; self.strain_ev(1) = self.init_strain;
           for i = 1:(length(self.t_list)-1)   
               self.residual_handle = @(ns) self.stress_ev(i+1) - ...
                   self.stressInc( cH, self.delta_t, ns, self.strain_ev(i) );
               self.strain_ev(i+1) = self.solve_at_time( self.t_list(i) ); 
               
               % As a check to see we don't have convergence for strain
               if isnan(self.strain_ev(i+1))
                    fprintf("\nTime %f :\tWARNING! NO CONVERGENCE! NEXT STRAIN IS NaN!", self.t_list(i));
               end
               
               cH = self.HInc(cH, self.delta_t, self.strain_ev(i+1), self.strain_ev(i));
               self.init_guess = self.strain_ev(i);
               fprintf('\nTime %f: Number of iterations = %d', self.t_list(i), self.num_iters);
               self.num_iters = 0;
               self.first_exec = true;
               self.first_residual = 9999;      % Reset first_residual value for use in next time step

               
           end

           fprintf('\n')
           
       end
       
       
       % Helper function to solve for strain
       function r = solve_at_time(self, t)
           % r is the returned value and we set it to the initial guess. 
           % increment is how much we increment/change r after each iteration
           r = self.init_guess;
           current_residual = 9999;
           % This computes a numerical derivative of the residual at the initial guess. If the modified
           % NR is used, then this is the derivative that will be used throughout solution finding.
           denom = (self.residual_handle(r+self.delta_x) -  ...
                    self.residual_handle(r-self.delta_x))/(2*self.delta_x);
           
           % We keep running this iteration scheme until the residual is
           % less than 0.025 AND the residual ratio is under 0.5. As we
           % approach a solution, the residual and residual ratio should
           % get smaller and smaller. The end criterion was kept numerically less
           % because otherwise we wouldn't get convergence for all
           % instances of delta_t. However this is countered by imposing
           % TWO end criteria that has to be met. It should be noted that
           % the strain computed is quite sensitive to bot the end 
           % criterion and the initial guess.
           
           while (abs(self.residual_handle(r))>self.abs_thresh) || ...
                 (abs(current_residual/self.first_residual) > self.rel_thresh)
             
                self.num_iters = self.num_iters + 1;
                if self.num_iters > 10000
                    r = NaN;
                    break;
                end

                numer = -self.residual_handle(r);    % this computes the residual at our guess

                if self.first_exec
                       self.first_residual = -numer;
                       self.first_exec = false;
                end        

                % This is to help log the residual at t = 0.5 to help us with the assignment
                % It also logs the respective number of iterations. This portion doesn't play a
                % role in finding the solution.
                if t == 0.5
                    self.residual_ratio = [self.residual_ratio, abs(numer/self.first_residual)];
                end

                % If we use the modified NR, then we don't recompute the derivative of the residual at
                % subsequent guesses. If we use the classic NR, then we recompute the derivative
                if ~self.modified
                    denom = (self.residual_handle(r+self.delta_x) -  ...
                            self.residual_handle(r-self.delta_x))/(2*self.delta_x);
                end
                
                increment = numer/denom;    % We increment our guess by -f/f', where f is the residual
                r = r + increment;          % Update our guess
                current_residual = self.residual_handle(r);
                
                
                if t == 0.5 && ~( (abs(self.residual_handle(r))>self.abs_thresh) || ...
                 (abs(current_residual/self.first_residual) > self.rel_thresh) )
                    self.residual_ratio = [self.residual_ratio, abs(current_residual/self.first_residual)];
             
                end

           end
           

                      
       end
       
   end    
    
end