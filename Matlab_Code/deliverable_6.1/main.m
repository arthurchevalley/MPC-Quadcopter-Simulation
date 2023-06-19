%--------------------------------------------------------------------------
%% ----------------------- MAIN of deliverable 6.1 ------------------------
% -------------------------------------------------------------------------
% Team members: - Pereira Portela Tifanny 
%               - Chevalley Arthur 
%               - Mermoud Paco
%
% Date: Autumn 2020
%
%
% The aim of this deliverable is to compute the non-linear MPC controller
% in order to follow a given path.
%
% -------------------------------------------------------------------------

clc

quad = Quad();

% Choice of non-linear MPC controller method
prompt = (['Which Non-linear method do you want to use ?',...
             '\nEnter 1 for the tifanny, 2 for the Arthur, anything else to stop\n']);
    NMPC_Choice = input(prompt);
    if NMPC_Choice == 1
        CTRL = My_ctrl_NMPC(quad);
    elseif  NMPC_Choice == 2 
        CTRL = ctrl_MPC(quad);
    else 
        error('Successful exit');
    end

% SImulate the non-linear MPC controller
sim = quad.sim(CTRL) 
quad.plot(sim)
