function [ Sxx ] = estimate_psd(X,lambda)

Sxx = filter(1-lambda, [1 -lambda], abs(X).^2, [], 2);