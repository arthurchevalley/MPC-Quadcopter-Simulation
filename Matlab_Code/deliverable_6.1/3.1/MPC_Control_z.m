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

      % SET THE HORIZON HERE
      N = 20;
      Q = 10*eye(n);
      R = 2*5*eye(m);
      
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
        M = [1;-1]; m = [0.3; 0.2];
        %state
        % x in X = { x | Fx <= f }
        % zero along z
        F=[];f=[];
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
%% Plot the sets  
        %{
        if(size(mpc.A,2)>2)
            figure('Name','Z subsystem invariant set','NumberTitle','off')
            subplot(1,3,1)
            Xf.projection(1:2).plot('alpha',0.9,'linewidth',1);
            hold on;
            subplot(1,3,2)
            Xf.projection(2:3).plot('alpha',0.9,'linewidth',1, 'colors','g');
            subplot(1,3,3)
            Xf.projection(3:4).plot('alpha',0.9,'linewidth',1, 'colors','b');
        else
            figure('Name','Z subsystem invariant set','NumberTitle','off')
            Xf.plot('alpha',0.2,'linewidth',3,'color','k')
        end
        %}
        [Ff,ff] = double(Xf);
        %% 4.1
      con = (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1)) + (M*u(:,1) <= m);
      obj = (u(:,1)-us)'*R*(u(:,1)-us);
    for i = 2:N-1
        con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i));
        con = con + (M*u(:,i) <= m);
        obj = obj + (x(:,i)-xs)'*Q*(x(:,i)-xs) + (u(:,i)-us)'*R*(u(:,i)-us);
    end
    con = con + (Ff*x(:,N) <= ff);
    obj = obj + (x(:,N)-xs)'*Qf*(x(:,N)-xs);

%% 5
      %{
        con = (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1) + mpc.B*d_est) + (M*u(:,1) <= m);
        obj = (u(:,1)-us)'*R*(u(:,1)-us);
        for i = 2:N-1
            con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i)+ mpc.B*d_est);
            con = con + (M*u(:,i) <= m);
            obj = obj + (x(:,i)-xs)'*Q*(x(:,i)-xs) + (u(:,i)-us)'*R*(u(:,i)-us);
        end
        con = con + (Ff*x(:,N) <= ff);
        obj = obj + (x(:,N)-xs)'*Qf*(x(:,N)-xs);
        %}
      
      
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
      % Constraints
      %input
      % u in U = { u | Hu <= k }
      Hu = [1;-1]; ku = [0.3; 0.2];
      %state
      % x in X = { x | Hx <= k }
      % no state constraint
      Hx = []; kx = [];
        
      Rs = eye(size(us,1));
      con = (xs == (mpc.A*xs + mpc.B*us)) + (ref == mpc.C*xs);      
      con = con + (Hu*us <= ku);
      obj = us'*Rs*us;
      

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
      
      %% 4.1
      
      [nx, nd] = size(mpc.B);
      A_bar = [];
      B_bar = [];
      C_bar = [];
      L = [];
      
      %% 5
      %{
      [nx, nd] = size(mpc.B);
      A_bar = [mpc.A mpc.B; zeros(1, nx) 1];
      B_bar = [mpc.B; zeros(nd)];
      C_bar = [mpc.C zeros(nd)];
      L = -place(A_bar', C_bar', [0.1, 0.11, 0.12])';
      %}
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
  end
end
