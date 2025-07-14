%GEODESIC EXAMPLE PLOT
% Geodesic paths on an n×n grid, single trial
clear; close all;

%% 0) Parameters
n = 20;    % grid size
t = 1;     % number of trials (we’ll only do one)

%% 1) Generate the random environment (one trial)
A = exprnd(1, [n, n]);

%% 2) Build DP tables for max and min passage times
dp = zeros(n, n);
dq = zeros(n, n);

% initialize (1,1)
dp(1,1) = A(1,1);
dq(1,1) = A(1,1);

% first row
for j = 2:n
    dp(1,j) = dp(1,j-1) + A(1,j);
    dq(1,j) = dq(1,j-1) + A(1,j);
end
% first column
for i = 2:n
    dp(i,1) = dp(i-1,1) + A(i,1);
    dq(i,1) = dq(i-1,1) + A(i,1);
end

% interior
for i = 2:n
    for j = 2:n
        dp(i,j) = max(dp(i-1,j), dp(i,j-1)) + A(i,j);
        dq(i,j) = min(dq(i-1,j), dq(i,j-1)) + A(i,j);
    end
end

%% 3) Backtrack the LAST-PASSAGE (max) geodesic from (n,n) to (1,1)
i = n; j = n;
path_dp = [i, j];
while i > 1 || j > 1
    if i > 1 && (j == 1 || dp(i-1,j) >= dp(i,j-1))
        i = i - 1;
    else
        j = j - 1;
    end
    path_dp(end+1, :) = [i, j];
end
path_dp = flipud(path_dp);

%% 4) Backtrack the FIRST-PASSAGE (min) geodesic
i = n; j = n;
path_dq = [i, j];
while i > 1 || j > 1
    if i > 1 && (j == 1 || dq(i-1,j) <= dq(i,j-1))
        i = i - 1;
    else
        j = j - 1;
    end
    path_dq(end+1, :) = [i, j];
end
path_dq = flipud(path_dq);

%% 5) Plot the two geodesics over the random field
figure;
imagesc(A);
colormap summer;
colorbar;
axis equal tight;
set(gca,'YDir','normal');  % so (1,1) is bottom-left
hold on;
% overlay paths (x = column j, y = row i)
plot(path_dp(:,2), path_dp(:,1), '-r', 'LineWidth', 2);
plot(path_dq(:,2), path_dq(:,1), '-b', 'LineWidth', 2);
legend('Last passage path', 'First passage path', 'Location', 'best');
title('Directed First and Last Passage Percolation on a Random Grid');
hold off;
