classdef MPC_Control_yaw < MPC_Control
  
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
      
      % Steady-state targets (Ignore this before Todo 3.2)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % -------------------PARAMS-------------------
      N = 5;
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

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      con = [];
      obj = 0;

      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      G = [1; -1];
      g = [0.2; 0.2];
     
      
      %-------------------- DROP TERMINAL SET ----------------------
      %{
      % Computing the infinite horizon LQR
      [Kih,P,~] = dlqr(mpc.A,mpc.B,Q_prob2,R_prob2,[]);
      Kih = -Kih; 

      % autonomous system A_star = A + B*K_ih
      A_cl = mpc.A + mpc.B*Kih;

      % Compute maximal invariant set
      V = [Kih; -Kih];

      Init = Polyhedron(V,g);
      
      Gamma = Init;
      
      compteur_yaw = 0;
      
      while 1 
          pre_Gamma = Polyhedron(double(Gamma.A)*A_cl, double(Gamma.b));
          intersection = Polyhedron([Gamma.A; pre_Gamma.A], [Gamma.b; pre_Gamma.b]);
          
          if intersection == Gamma
              break 
          end
          compteur_yaw = compteur_yaw + 1
          Gamma = intersection;   
      end 
      
      terminal_set = Gamma;
      %}
      
      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      con = [];
      obj = 0;
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % prob 1
      for i = 1:N-1
        con = [con, x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)]; % System dynamics
        con = [con, G*u(:,i) <= g]; % Input constraints
        obj = obj + (x(:,i)-xs)'*Q_prob1*(x(:,i)-xs) + (u(:,i)-us)'*R_prob1*(u(:,i)-us); % Cost function
      end
      
      %-------------------- DROP TERMINAL SET ----------------------
      % con = [con, terminal_set.A*(x(:,N)-xs) <= terminal_set.b]; % Terminal constraint
      % con = [con, F*(x(:,N)-xs) <= f]; % Terminal constraint
     
      % obj = obj + (x(:,N)-xs)'*P*(x(:,N)-xs); % Terminal weight
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
      
      % Reference position (Ignore this before Todo 3.2)
      ref = sdpvar;            
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE       
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      con = [];
      obj = 0;

      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      G = [1; -1];
      g = [0.2; 0.2];
      
      con = [G*us <= g, xs == mpc.A*xs + mpc.B*us, ref == mpc.C*xs + mpc.D];
      obj = us^2;
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
