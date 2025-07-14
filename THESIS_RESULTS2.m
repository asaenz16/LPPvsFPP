%PQD SIMULATION
clear;

n = 100; % Grid size
t = 10000; % Number of samples


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

m = 300;  
x_vals = linspace(min(dp_values), max(dp_values), m);
y_vals = linspace(min(dq_values), max(dq_values), m);

Fxy = zeros(m,m);   % joint CDF
Fx  = zeros(1,m);   % marginal CDF of dp
Fy  = zeros(1,m);   % marginal CDF of dq

% 3) fill in Fx and Fy
for i = 1:m
    Fx(i) = mean( dp_values <= x_vals(i) );
    Fy(i) = mean( dq_values <= y_vals(i) );
end

% 4) fill in joint F(x,y)
for i = 1:m
  for j = 1:m
    Fxy(j,i) = mean( dp_values <= x_vals(i)  &  dq_values <= y_vals(j) );
  end
end

Delta = Fxy - (Fy' * Fx);

figure
surf(x_vals, y_vals, Delta)
xlabel('x'); ylabel('y'); zlabel('\Delta F = F(x,y)-F(x)F(y)')
title('Empirical Cumulative Covariance Function')
shading interp
colorbar