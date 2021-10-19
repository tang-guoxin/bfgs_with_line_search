function lambda = linear_search(fun, gfun, x0, alpha, beta)
%% ���ݷ��������õ����Ų���
if nargin == 3
   alpha = rand * 0.5;
   beta = rand;
end
%% 
t = 1;
dert_x = -gfun(x0);
% fv = fun(x0 + t*dert_x);
while fun(x0 + t * dert_x) > fun(x0) + alpha * t * gfun(x0)' * dert_x
   t = beta * t; 
end
%% 
lambda = t;
end
