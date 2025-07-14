%COVARIANCE CALCULATIONS
clear;

n = 2; % Grid size
t = 1000000; % Number of samples


% Storage for passage times
dp_values = zeros(1, t); % For storing dp (maximal passage times)
dq_values = zeros(1, t); % For storing dp (maximal passage times)
for k = 1:t
    % Generate random environment in one step
    A = exprnd(1, [n, n]); 
    
    % Compute `dp` for MAXIMAL path weights
    dp = zeros(n, n);
    dp(1, :) = cumsum(A(1, :)); % First row cumulative sum
    dp(:, 1) = cumsum(A(:, 1)); % First column cumulative sum
    % Compute `dq` for MINIMAL path weights
    dq = zeros(n, n);
    dq(1, :) = cumsum(A(1, :)); % First row cumulative sum
    dq(:, 1) = cumsum(A(:, 1)); % First column cumulative sum
    for i = 2:n
        for j = 2:n
            dp(i, j) = max(dp(i, j-1), dp(i-1, j)) + A(i, j);
            dq(i, j) = min(dq(i, j-1), dq(i-1, j)) + A(i, j);
        end
    end
    
    % Store the final passage times
    dp_values(k) = dp(n, n);
    dq_values(k) = dq(n, n);
end


% === Compute Moments ===
max_moment = 1;
moments_dp = zeros(1, max_moment);
moments_dq = zeros(1, max_moment);
moments_product = zeros(1, max_moment);

for k = 1:max_moment
    moments_dp(k) = mean(dp_values.^k);
    moments_dq(k) = mean(dq_values.^k);
    moments_product(k) = mean((dp_values .* dq_values).^k);
end
Indep_test=moments_product-moments_dp.*moments_dq;



