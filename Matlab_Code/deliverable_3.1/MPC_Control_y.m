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
      
      % Steady-state targets (Ignore this before Todo 3.2)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % -------------------PARAMS-------------------
      N = 20; % horizon length 
      
      % parameters of the finite horizon problem
      Q_prob1 = 10*eye(4);
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
     
      % Computing the infinite horizon LQR Gain
      [Kih,P,~] = dlqr(mpc.A,mpc.B,Q_prob2,R_prob2,[]);
      Kih = -Kih; 

      % Closed-loop system 
      A_cl = mpc.A + mpc.B*Kih;

      % Compute maximal invariant set
      V = [Kih; -Kih];

      Gamma = Polyhedron([F;V], [f;g]);
      
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
      subplot(1,3,1) 
      terminal_set.projection(1:2).plot('alpha', 0.1)
      title('Terminal set projection on dim 1:2')
      xlabel('Roll angular velocity in [rad/s]')
      ylabel('Roll angle in [rad]')
      
      subplot(1,3,2) 
      terminal_set.projection(2:3).plot('alpha', 0.1)
      title('Terminal set projection on dim 2:3')
      xlabel('Roll angle in [rad]')
      ylabel('linear velocity along y in [m/s]')
      
      subplot(1,3,3)
      terminal_set.projection(3:4).plot('alpha', 0.1)
      title('Terminal set projection on dim 3:4')
      xlabel('linear velocity along y in [m/s]')
      ylabel('position along y in [m]')
      
      sgt = sgtitle('Projections of terminal set for system y','Color','black');
      sgt.FontSize = 20;

      % Constraints and objective of the finite horizon problem
      con = [];
      obj = 0;
     
      for i = 1:N-1
        con = [con, x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)]; % System dynamics
        con = [con, F*x(:,i) <= f]; % State constraints
        con = [con, G*u(:,i) <= g]; % Input constraints
        obj = obj + x(:,i)'*Q_prob1*x(:,i) + u(:,i)'*R_prob1*u(:,i); % Cost function
      end
      
      con = [con, terminal_set.A*x(:,N) <= terminal_set.b]; % Terminal constraint
      obj = obj + x(:,N)'*P*x(:,N); % Terminal weight 
      
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
      
      % Reference position, used in a future deliverable
      ref = sdpvar;            
            
      % Constraints and objectives, done in a future deliverable
      con = [];
      obj = 0;
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
