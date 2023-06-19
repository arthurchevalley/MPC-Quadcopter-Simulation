classdef MPC_Control_z < MPC_Control
  properties
    A_bar, B_bar, C_bar % Augmented system for disturbance rejection    
    L                   % Estimator gain for disturbance rejection
  end
  
  methods
    function mpc = MPC_Control_z(sys, Ts)
      mpc = mpc@MPC_Control(sys, Ts);
      
      [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
    end
    
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   d_est  - disturbance estimate
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % Disturbance estimate
      d_est = sdpvar(1);

      % -------------------PARAMS-------------------
      
      N = 10; % Horizon length
      
      % parameters of the finite horizon problem
      Q_prob1 = 10*eye(2);
      R_prob1 = 1;
      
      % parameters of the infinite horizon problem
      Q_prob2 = Q_prob1;
      R_prob2 = R_prob1;
      % ---------------------------------------------
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      ny = size(mpc.C,1);
      
      % Input constraints
      G = [1; -1];
      g = [0.3; 0.2];
      
      % Matrix handling the bias
      Bd = mpc.B;

      
      % Constraints and objective of the finite horizon problem
      con = [];
      obj = 0;
      
      for i = 1:N-1
        con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i) + Bd*d_est); % System dynamics
        con = con + (G*u(:,i) <= g); % Input constraints
        obj = obj + (x(:,i)-xs)'*Q_prob1*(x(:,i)-xs) + (u(:,i)-us)'*R_prob1*(u(:,i)-us); % Cost function
      end
      
      obj = obj + (x(:,N)-xs)'*Q_prob1*(x(:,N)-xs); % Terminal weight 

      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us, d_est}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      %   d_est  - disturbance estimate
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      ny = size(mpc.C,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position 
      ref = sdpvar;
            
      % Disturbance estimate 
      d_est = sdpvar(1);
      
      % Input constraints
      G = [1; -1];
      g = [0.3; 0.2];
      
      % Bias matrices
      Bd = mpc.B;
      Cd = zeros(ny,1);
      
      % Constraints and objectives of the steady-state subject to the bias
      con = [G*us <= g, xs == mpc.A*xs + mpc.B*us + Bd*d_est, ref == mpc.C*xs + Cd*d_est];
      obj = us^2;

      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
    end
    
    
    % Compute augmented system and estimator gain for input disturbance rejection
    function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
      
      %%% Design the matrices A_bar, B_bar, L, and C_bar
      %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
      %%% converges to the correct state and constant input disturbance
      %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      % check if the system is observable 
      
      nx = size(mpc.A,1);
      nu = size(mpc.B,2);
      ny = size(mpc.C,1);
     
      % Bias matrices
      Bd = mpc.B;
      Cd = zeros(ny,1);

      % Check if the augmented system is observable 
      Q = [mpc.C;mpc.C*mpc.A];
      test_rank = [mpc.A-eye(2), Bd; mpc.C, Cd];
      if (rank(Q) == nx) & (rank(test_rank) == (nu + nx))
          disp('The augmented system is observable');
      end 

      % Observer matrices 
      A_bar = [mpc.A, Bd; zeros(1,nx), 1];
      B_bar = [mpc.B; zeros(1, nu)];
      C_bar = [mpc.C, Cd];
      % Observer gain
      L = -place(A_bar',C_bar',[0.4,0.5,0.6])';

    end

    
  end
end
