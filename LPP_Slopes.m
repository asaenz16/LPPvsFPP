clear;
n = 100;
p = 0.5;
t = 5000;
F_value = 30;
L_value = 180;
data = sparse(2*n, 2*n); % Preallocate as sparse

for k = 1:t
    % Generate random environment in one step
    A = binornd(1, p, [n, n]); 
    
    % Compute dp for MAXIMAL path weights
    dp = zeros(n, n);
    dp(1, :) = cumsum(A(1, :)); % First row cumulative sum
    dp(:, 1) = cumsum(A(:, 1)); % First column cumulative sum
    
    for i = 2:n
        for j = 2:n
            dp(i, j) = max(dp(i, j-1), dp(i-1, j)) + A(i, j);
        end
    end
    
    % Compute dq for MINIMAL path weights
    dq = zeros(n, n);
    dq(1, :) = cumsum(A(1, :)); % First row cumulative sum
    dq(:, 1) = cumsum(A(:, 1)); % First column cumulative sum
    
    for i = 2:n
        for j = 2:n
            dq(i, j) = min(dq(i, j-1), dq(i-1, j)) + A(i, j);
        end
    end
    
    % Update data matrix
    dq_index = max(1, dq(n, n));
    dp_index = max(1, dp(n, n));
    data(dq_index, dp_index) = data(dq_index, dp_index) + 1;
end

% Convert sparse matrix to full for plotting
data = full(data);
data = data / t; % Normalize to a pmf: P(L and F)
marginal_dataL = sum(data); % P(L)
marginal_dataF = sum(data'); % P(F)
conditional_dataL = zeros(size(data)); % P(L|F)
for i = 1:2*n
    if marginal_dataF(i) ~= 0
        conditional_dataL(i, :) = data(i, :) ./ marginal_dataF(i);
    end
end
conditional_dataF = zeros(size(data)); % P(F|L)
for i = 1:2*n
    if marginal_dataL(i) ~= 0
        conditional_dataF(:, i) = data(:, i) ./ marginal_dataL(i);
    end
end

% Define the range for rows and columns
rowRange = 1:60;
colRange = 140:2*n;

% Extract the submatrices for the selected range
subDataL = conditional_dataL(rowRange, colRange);
subDataF = conditional_dataF(rowRange, colRange);

% Create the subplots for specific densities
figure;

% --------------------------
% Plot 1: P(L | F=30)
% --------------------------
subplot(2, 1, 1); % 2 rows, 1 column, first plot

imagesc(colRange, rowRange, subDataL); % Use colRange and rowRange for axis labels
colorbar; % Add a color bar to indicate scale
xlabel('L values');
ylabel('F values');
title('Heatmap of P(L | F)');
set(gca, 'YDir', 'normal'); % Ensure correct orientation


subplot(2, 1, 2); % 2 rows, 1 column, second plot
subplot(2, 1, 2); % 2 rows, 1 column, second plot
imagesc(colRange, rowRange, subDataF); % Use colRange and rowRange for axis labels
colorbar; % Add a color bar to indicate scale
xlabel('L values');
ylabel('F values');
title('Heatmap of P(F | L)');
set(gca, 'YDir', 'normal'); % Ensure correct orientation
