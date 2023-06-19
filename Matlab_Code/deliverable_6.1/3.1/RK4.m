function [x_next] = RK4(X,U,h,quad)
%
% Inputs : 
%    X, U current state and input
%    h    sample period
%    f    continuous time dynamics f(x,u)
% Returns
%    State h seconds in the future
%

% Runge-Kutta 4 integration
% write your function here
   k1 = quad.f(X,U);
   k2 = quad.f(X+h/2*k1,U);
   k3 = quad.f(X+h/2*k2,U);
   k4 = quad.f(X+h*k3,U);
   x_next = X + h/6*(k1+2*k2+2*k3+k4);
