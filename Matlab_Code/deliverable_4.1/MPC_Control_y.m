classdef MPC_Control_y < MPC_Control
  
  methods
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % -------------------PARAMS-------------------
      N = 14; % Horizon length
      
      % parameters of the finite horizon problem
      Q_prob1 = 100*eye(4);
      R_prob1 = 1;
      
      % parameters of the infinite horizon problem
      Q_prob2 = Q_prob1;
      R_prob2 = R_prob1;
      % ---------------------------------------------
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      % State and input constraints 
      F = [0 1 0 0 ; 0 -1 0 0];
      f = [0.035 ; 0.035];
      G = [1; -1];
      g = [0.3; 0.3];
        
    
      % Constraints and objective of the finite horizon problem
      con = [];
      obj = 0;

      for i = 1:N-1
        con = [con, x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)]; % System dynamics
        con = [con, F*x(:,i) <= f]; % State constraints
        con = [con, G*u(:,i) <= g]; % Input constraints
        obj = obj + (x(:,i)-xs)'*Q_prob1*(x(:,i)-xs) + (u(:,i)-us)'*R_prob1*(u(:,i)-us); % Cost function
      end
      
      con = [con, F*x(:,N) <= f]; % Terminal constraint
      obj = obj + (x(:,N)-xs)'*Q_prob1*(x(:,N)-xs); % Terminal weight
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position 
      ref = sdpvar;   
      
      % State and input constraints 
      F = [0 1 0 0 ; 0 -1 0 0];
      f = [0.035 ; 0.035];
      G = [1; -1];
      g = [0.3; 0.3];
      
      % Constraints and objective for the steady-state computation
      con = [G*us <= g, F*xs <= f, xs == mpc.A*xs + mpc.B*us, ref == mpc.C*xs + mpc.D];
      obj = us^2;
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
