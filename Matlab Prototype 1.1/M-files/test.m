
Q = [1 0 0 0; 1 1 0 0; 2 1 2 0; 1 1 1 1];
P = Q*Q.'; 

[M,H] = optimize_trace(P,2);

tr1 = trace(M); 

tr2 = trace(P1) - 2*h.'*h;

