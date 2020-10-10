function [x_j, error] = Jacobi(A,x0,b, tolerance)
    error = tolerance + 1;
    M = -triu(A, 1);
    N = -tril(A, -1);
    D = diag(diag(A));
    D_inv = inv(D);
    T = D_inv*(M + N);
    c = D_inv*b;
    x_i = x0;
    while error>tolerance
      x_j = T*x_i + c;
      error = norm(x_j - x_i, Inf);  
      x_i = x_j;
    endwhile
    
endfunction