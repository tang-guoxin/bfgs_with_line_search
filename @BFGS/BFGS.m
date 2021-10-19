classdef BFGS
    %% BFGS 算法
    %% 公共参数
    properties
        fval = inf;
        eps = 1e-6;
        opt = 'BFGS';
        x0
        fun
        gfun
        step
        max_iter
        dim
        limit
        alpha
        beta
    end
    %% 方法
    methods
        function this = BFGS(fun, gfun, x0, options)
           %% 初始化函数
            this.fun = fun;
            this.gfun = gfun;
            this.x0 = x0;
            this.dim = length(x0);
            if nargin == 3
                this.step = 'linear-search';
                this.max_iter = 100;
                this.opt = 'BFGS';
                this.limit = 10;
                this.alpha = -1;
                this.beta = -1;
            else
                this.step = options.step;
                this.max_iter = options.max_iter;
                this.opt = options.opt;
                this.limit = options.limit;
                this.alpha = options.alpha;
                this.beta = options.beta;
            end
        end
       %% 
        function results = fit(this, verbose)
           %% 这个是主函数, 默认不显示细节, 设置verbose=true来显示迭代过程
            if nargin == 1
                verbose = false;
            end
            tic;
            kcount = 0;
            Bk = eye(this.dim);
            gk = this.gfun(this.x0);
            % 内存限制矩阵
            sk_mat = zeros(this.limit, this.dim);
            yk_mat = zeros(this.limit, this.dim);
            %  主循环
            for k = 1:this.max_iter
                if strcmp('L-BFGS', this.opt)
                    % dk = limit_memory(sk, yk) * gk';
                    dk = -inv(Bk) * gk';
                else
                    dk = -inv(Bk) * gk';
                end
                if strcmp(this.step, 'linear-search')
                    if this.alpha == -1 || this.beta == -1
                        lambda = linear_search(this.fun, this.gfun, this.x0);
                    else
                        lambda = linear_search(this.fun, this.gfun, this.x0, this.alpha, this.beta);
                    end
                else
                    lambda = this.step;
                end
                sk = lambda * dk;
                this.x0 = this.x0 + sk;
                gk1 = this.gfun(this.x0)';
                if norm(gk1) < this.eps
                    break;
                end
                kcount = kcount + 1;
                yk = gk1 - gk;
                gk = gk1;
                if strcmp('L-BFGS', this.opt)
                    if kcount <= this.limit
                        sk_mat(kcount, :) = sk;
                        yk_mat(kcount, :) = yk;
                    else
                        sk_mat(1:end-1, :) = sk_mat(2:end, :);
                        sk_mat(end, :) = sk;
                        yk_mat(1:end-1, :) = yk_mat(2:end, :);
                        yk_mat(end, :) = yk;
                    end
                else
                   Bk = Bk + (yk * yk') / (yk' * sk) - (Bk * sk *  sk * Bk) / (sk' * Bk * sk);
                end
                if verbose == true
                    fprintf('iter = %d \t x = %f \t fval = %f \t gradf = %f \n', k, this.x0, this.fval, this.gfun(this.x0))
                end
                this.fval = this.fun(this.x0);
            end
           %% return
            results.x0 = this.x0;
            results.fval = this.fval;
            results.iteration = kcount;
            results.cost_time = toc;
        end
    end
end


function lambda = linear_search(fun, gfun, x0, alpha, beta)
%% 回溯法线搜索得到最优步长
if nargin == 3
   alpha = rand * 0.5;
   beta = rand;
end
%% 
t = 1;
dert_x = -gfun(x0);
while fun(x0 + t * dert_x) > fun(x0) + alpha * t * gfun(x0)' * dert_x
   t = beta * t; 
end
%% 
lambda = t;
end