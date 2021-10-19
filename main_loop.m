tic, clc, clear, close all, format longG
%% 
f = @(x) x.^2 + 5*sin(x);
gf = @(x) 2*x + 5*cos(x);
%% 
subplot 211
ezplot(f)
legend({'ԭ����'})
subplot 212
ezplot(gf)
legend({'һ�׵���'})

%% 
options.opt = 'BFGS';               % �Ż�����, Ĭ��BFGS, L-BFGS��û�����
options.step = 'linear-search';     % �Ƿ����������
options.max_iter = 100;             % ����������
options.limit = 10;                 % �ڴ����Ʊ������ʷ��������
options.alpha = -1;                 % ���������� alpha, ����Ϊ -1 �����ý�ʹ��������� alpha \in (0, 0.5)
options.beta = -1;                  % ���������� beta, ����Ϊ -1 �����ý�ʹ��������� beta \in (0, 1)
%% ���� BFGS �㷨
model = BFGS(f, gf, 2, options)
% ��ʼִ���㷨 ����Ϊ false ʱ�����������������, ����һ���ṹ��, �������Ž�, ����ֵ, �������������ʱ��
res = model.fit(true)

% res.x0
% res.fval

%%
toc
