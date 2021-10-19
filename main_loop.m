tic, clc, clear, close all, format longG
%% 
f = @(x) x.^2 + 5*sin(x);
gf = @(x) 2*x + 5*cos(x);
%% 
subplot 211
ezplot(f)
legend({'原函数'})
subplot 212
ezplot(gf)
legend({'一阶导数'})

%% 
options.opt = 'BFGS';               % 优化方法, 默认BFGS, L-BFGS还没有完成
options.step = 'linear-search';     % 是否进行线搜索
options.max_iter = 100;             % 最大迭代次数
options.limit = 10;                 % 内存限制保存的历史向量个数
options.alpha = -1;                 % 线搜索参数 alpha, 设置为 -1 或不设置将使用随机参数 alpha \in (0, 0.5)
options.beta = -1;                  % 线搜索参数 beta, 设置为 -1 或不设置将使用随机参数 beta \in (0, 1)
%% 创建 BFGS 算法
model = BFGS(f, gf, 2, options)
% 开始执行算法 参数为 false 时将不再输出迭代过程, 返回一个结构体, 包括最优解, 最优值, 迭代次数与迭代时间
res = model.fit(true)

% res.x0
% res.fval

%%
toc
