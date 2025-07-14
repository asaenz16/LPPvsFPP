%GEODESIC SURF
clear;

%% 0) Parameters
n = 100;     % number of columns
t = 1000000;    % number of trials

%% 1) Preallocate storage for the drop-columns
dp_jump = zeros(1, t);   % last-passage (max) path
dq_jump = zeros(1, t);   % first-passage (min) path

for k = 1:t
    %% 2) Generate 2×n environment
    A = exprnd(1, [2, n]);

    %% 3) Build 2×n DP tables
    dp = zeros(2, n);
    dq = zeros(2, n);
    dp(1, :) = cumsum(A(1, :));
    dp(:, 1) = cumsum(A(:, 1));
    dq(1, :) = cumsum(A(1, :));
    dq(:, 1) = cumsum(A(:, 1));
    for j = 2:n
        dp(2, j) = max(dp(1, j), dp(2, j-1)) + A(2, j);
        dq(2, j) = min(dq(1, j), dq(2, j-1)) + A(2, j);
    end

    %% 4) Backtrack the MAXIMAL geodesic from (2,n) → (1,1)
    i = 2; j = n;
    path = [i, j];
    while i > 1 || j > 1
        if i > 1 && (j == 1 || dp(i-1, j) >= dp(i, j-1))
            i = i - 1;
        else
            j = j - 1;
        end
        path(end+1, :) = [i, j];
    end
    path = flipud(path);  % now from (1,1) → (2,n)
    first2 = find(path(:,1) == 2, 1, 'first');
    dp_jump(k) = path(first2, 2);

    %% 5) Backtrack the MINIMAL geodesic from (2,n) → (1,1)
    i = 2; j = n;
    path = [i, j];
    while i > 1 || j > 1
        if i > 1 && (j == 1 || dq(i-1, j) <= dq(i, j-1))
            i = i - 1;
        else
            j = j - 1;
        end
        path(end+1, :) = [i, j];
    end
    path = flipud(path);
    first2 = find(path(:,1) == 2, 1, 'first');
    dq_jump(k) = path(first2, 2);
end

%% 6) Empirical joint PMF over {1,…,n}×{1,…,n}
% Build a two-column subscript array, then count
subs = [dq_jump.' , dp_jump.'];    % rows=dq_jump, cols=dp_jump
counts = accumarray(subs, 1, [n, n]);
jointPMF = counts / t;   % sums to 1

%% 7) Plot (normalized to [0,1]^2 with density scaling)
[xGrid, yGrid] = meshgrid((1:n)/n, (1:n)/n);  % normalized domain
jointPDF = jointPMF * n^2;                   % scale to density

figure
surf(xGrid, yGrid, jointPDF)
xlabel('Last Passage Geodesic (k/N)')
ylabel('First Passage Geodesic (j/N)')
zlabel('P(k/N=x, j/N=y)')
title('Joint PMF of Geodesic Position on 2×N Lattice')
shading interp
colorbar

