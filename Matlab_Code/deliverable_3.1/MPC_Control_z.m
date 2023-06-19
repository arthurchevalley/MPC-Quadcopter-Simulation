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
      
      % Steady-state targets, used in a future deliverable
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % Disturbance estimate, used in a future deliverable
      d_est = sdpvar(1);

      % -------------------PARAMS-------------------
      
      N = 10; % horizon length
      
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
      
      % Input constraints
      G = [1; -1];
      g = [0.3; 0.2];
     
      % Computing the infinite horizon LQR Gain
      [Kih,P,~] = dlqr(mpc.A,mpc.B,Q_prob2,R_prob2,[]);
      Kih = -Kih; 

      % Closed-loop system
      A_cl = mpc.A + mpc.B*Kih;

      % Compute maximal invariant set
      V = [Kih; -Kih];

      Gamma = Polyhedron(V,g);
      
      while 1 
          pre_Gamma = Polyhedron(double(Gamma.A)*A_cl, double(Gamma.b));
          intersection = Polyhedron([Gamma.A; pre_Gamma.A], [Gamma.b; pre_Gamma.b]);
          
          if intersection == Gamma
              break 
          end
          Gamma = intersection;   
      end 
      terminal_set = Gamma;
     
      % Plot the terminal set 
      figure()
      terminal_set.plot('alpha', 0.1)
      xlabel('linear velocity along z in [m/s]')
      ylabel('position along z in [m]')
      sgt = sgtitle('Terminal set for system z','Color','black');
      sgt.FontSize = 20;

      % Constraints and objective of the finite horizon problem
      con = [];
      obj = 0;

      for i = 1:N-1
        con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)); % System dynamics
        con = con + (G*u(:,i) <= g); % Input constraints
        obj = obj + x(:,i)'*Q_prob1*x(:,i) + u(:,i)'*R_prob1*u(:,i); % Cost function
      end
      
      con = con + (double(terminal_set.A)*x(:,N) <= double(terminal_set.b)); % Terminal constraint
      obj = obj + x(:,N)'*P*x(:,N); % Terminal weight 

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
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position, used in a future deliverable
      ref = sdpvar;
            
      % Disturbance estimate, used in a future deliverable
      d_est = sdpvar(1);
      
      % Constraints and objectives, done in a future deliverable
      con = [];
      obj = 0;
      
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
    end
    
    
    % Compute augmented system and estimator gain for input disturbance rejection
    function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
      
      %%% Design the matrices A_bar, B_bar, L, and C_bar
      %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
      %%% converges to the correct state and constant input disturbance
      %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
      
      % Observer matrices, done in a future deliverable
      A_bar = [];
      B_bar = [];
      C_bar = [];
      L = [];
    end

    
  end
end
