function dk = limit_memory(sk_mat, yk_mat, kcount)
[m, ~] = size(sk_mat);
%% �����ڴ�: �ݹ� ?
ite = min(kcount, m);
for k = 1:ite
    sk = sk_mat(k, :);
    yk = yk_mat(k, :);
    %
    break
end
%% 
dk = 0;
end

