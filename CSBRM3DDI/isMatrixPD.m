function PD = isMatrixPD(M)
    [~,flag] = chol(M);  % if flag = 0, the input matrix is symmetric positive definite and the factorization was successful
    PD = ~ flag;
end