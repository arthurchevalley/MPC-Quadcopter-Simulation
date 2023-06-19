clc

quad = Quad();
CTRL = ctrl_MPC(quad);
sim = quad.sim(CTRL);
quad.plot(sim)