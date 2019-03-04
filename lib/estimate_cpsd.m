function [ Sx1x2 ] = estimate_cpsd(X1,X2,lambda)

Sx1x2 = filter(1-lambda, [1 -lambda], X1.*conj(X2), [], 2);