%COVARIANCE CALCULATIONS
clear;
meansvars=zeros([2,15]);
for idx=1:15
n = 100*idx; % Grid size
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


% === Compute mean and variance of dq_values ===
meansvars(1,idx) = mean(dq_values);
meansvars(2,idx)  = var(dq_values);
end

N = (1:15)*100;

figure;

% Left plot: meanvars(1,:)./N vs N
subplot(1,2,1);
plot(N, meansvars(1,:)./N, 'o-','LineWidth',1.5);
xlabel('N');
ylabel('Mean / N');
title('Normalized Mean');
grid on;

% Right plot: log(vars)/log(N) vs N
subplot(1,2,2);
plot(N, log(meansvars(2,:))./log(N), 'o-','LineWidth',1.5);
xlabel('N');
ylabel('log(Var) / log(N)');
title('Empirical Scaling Exponent');
grid on;


set([subplot(1,2,1), subplot(1,2,2)], 'XLim', [0 1500], 'YLim', [0 1]);
h = [subplot(1,2,1), subplot(1,2,2)];
set(h, 'XLim', [0 1500], 'YLim', [0 1]);

% add a blue dashed line at y = 0.6666667 on the 2nd subplot,
% give it a DisplayName, then show the legend
yl = yline(h(2), 0.6666667, 'b--', 'DisplayName','Conjecture', 'y = 0.6667');
legend(h(2), 'show');

