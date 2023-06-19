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
      
      % SET THE HORIZON HERE
      N = 20;
      Q = eye(n);
      R = 10*eye(m);
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      % Constraints
        %input
        % u in U = { u | Mu <= m }
        M = [1;-1]; m = [0.2; 0.2];
        %state
        % x in X = { x | Fx <= f }
        % zero for yaw
        F = [];f=[];
        % Compute LQR controller for unconstrained system
        [K,Qf,~] = dlqr(mpc.A,mpc.B,Q,R);
        % MATLAB defines K as -K, so invert its signal
        K = -K; 

        % Compute maximal invariant set
        Xf = polytope([F;M*K],[f;m]);
        Acl = [mpc.A+mpc.B*K];
        while 1
            prevXf = Xf;
            [T,t] = double(Xf);
            preXf = polytope(T*Acl,t);
            Xf = intersect(Xf, preXf);
            if isequal(prevXf, Xf)
                break
            end
        end
        
%% Plot the set
%{
        if(size(mpc.A,2)>2)
            figure('Name','Yaw subsystem invariant set','NumberTitle','off')
            subplot(1,3,1)
            Xf.projection(1:2).plot('alpha',0.9,'linewidth',1);
            hold on;
            subplot(1,3,2)
            Xf.projection(2:3).plot('alpha',0.9,'linewidth',1, 'colors','g');
            subplot(1,3,3)
            Xf.projection(3:4).plot('alpha',0.9,'linewidth',1, 'colors','b');
        else
            figure('Name','Yaw subsystem invariant set','NumberTitle','off')
            Xf.plot('alpha',0.2,'linewidth',3,'color','k')
        end
 %}       
        [Ff,ff] = double(Xf);
        
      con = (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1)) + (M*u(:,1) <= m);
      obj = u(:,1)'*R*u(:,1);
    for i = 2:N-1
        con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i));
        con = con + (M*u(:,i) <= m);
        obj = obj + (x(:,i)-xs)'*Q*(x(:,i)-xs) + (u(:,i)-us)'*R*(u(:,i)-us);
    end
    con = con + (Ff*x(:,N) <= ff);
    obj = obj + (x(:,N)-xs)'*Qf*(x(:,N)-xs);

      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
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
      % Constraints
      %input
      % u in U = { u | Mu <= m }
      Hu = [1;-1]; ku = [0.2; 0.2];
      %state
      % x in X = { x | Fx <= f }
      % for beta <= 0.035
      %Hx = []; kx = [];
        
      Rs = eye(size(us,1));
      con = (xs == (mpc.A*xs + mpc.B*us)) + (ref == mpc.C*xs);      
      con = con + (Hu*us <= ku);
      obj = us'*Rs*us;

      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
