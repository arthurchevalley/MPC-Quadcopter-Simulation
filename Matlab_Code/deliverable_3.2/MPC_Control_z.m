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
      
      % Steady-state targets (Ignore this before Todo 3.3)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);

      % -------------------PARAMS-------------------
      
      N = 5;
      
      % avant deli 41
      % N = 3
      % Q_prob1 = 10*eye(2);
      % R_prob1 = 1;
      
      Q_prob1 = 10*eye(2);
      R_prob1 = 1;
      
      Q_prob2 = Q_prob1;
      R_prob2 = R_prob1;
      % ---------------------------------------------
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      G = [1; -1];
      g = [0.3; 0.2];
      
      %-------------------- DROP TERMINAL SET ----------------------
      %{
     
      % Computing the infinite horizon LQR
      [Kih,P,~] = dlqr(mpc.A,mpc.B,Q_prob2,R_prob2,[]);
      Kih = -Kih; 

      % autonomous system A_star = A + B*K_ih
      A_cl = mpc.A + mpc.B*Kih;

      % Compute maximal invariant set
      % V = G * K_ih
      V = [Kih; -Kih];

      Init = Polyhedron(V,g);
      
      Gamma = Init;
      compteur_z = 0;
      
      while 1 
          pre_Gamma = Polyhedron(double(Gamma.A)*A_cl, double(Gamma.b));
          intersection = Polyhedron([Gamma.A; pre_Gamma.A], [Gamma.b; pre_Gamma.b]);
          
          if intersection == Gamma
              break 
          end
          compteur_z = compteur_z + 1
          Gamma = intersection;   
      end 
      
      terminal_set = Gamma;
     
      %}
      
      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      con = [];
      obj = 0;
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % pas de cost pour être en x0?
      % obj = (u(:,1)-us)'*R_prob1*(u(:,1)-us);
      % con = (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1)); % System dynamics
      % con = con + (G*u(:,1) <= g); % Input constraints
      
      % prob 1
      for i = 1:N-1
        con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)); % System dynamics
        con = con + (G*u(:,i) <= g); % Input constraints
        obj = obj + (x(:,i)-xs)'*Q_prob1*(x(:,i)-xs) + (u(:,i)-us)'*R_prob1*(u(:,i)-us); % Cost function
      end
      
      %-------------------- DROP TERMINAL SET ----------------------
      % con = [con, terminal_set.A*(x(:,N)-xs) <= terminal_set.b]; % Terminal constraint
      % con = [con, F*(x(:,N)-xs) <= f]; % Terminal constraint
      
      % obj = obj + (x(:,N)-xs)'*P*(x(:,N)-xs); % Terminal weight
      obj = obj + (x(:,N)-xs)'*Q_prob1*(x(:,N)-xs); % Terminal weight

      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
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
      
      % Reference position (Ignore this before Todo 3.3)
      ref = sdpvar;
            
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      G = [1; -1];
      g = [0.3; 0.2];
      
      con = [G*us <= g, xs == mpc.A*xs + mpc.B*us, ref == mpc.C*xs + mpc.D];
      obj = us^2;
      
      %con = [G*us <= g, xs == mpc.A*xs + mpc.B*us];
      %obj = (mpc.C*xs - ref)^2;

      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
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
 
      A_bar = [];
      B_bar = [];
      C_bar = [];
      L = [];
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
  end
end
