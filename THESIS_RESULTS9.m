clear; close all;
n = 20;
t = 2;

figure('Units','normalized','Position',[0.1 0.1 0.8 0.5]);
for k = 1:t
    A = exprnd(1, [n, n]);
    dp = zeros(n, n);
    dq = zeros(n, n);
    dp(1,1) = A(1,1);
    dq(1,1) = A(1,1);
    for j = 2:n
        dp(1,j) = dp(1,j-1) + A(1,j);
        dq(1,j) = dq(1,j-1) + A(1,j);
    end
    for i = 2:n
        dp(i,1) = dp(i-1,1) + A(i,1);
        dq(i,1) = dq(i-1,1) + A(i,1);
    end
    for i = 2:n
        for j = 2:n
            dp(i,j) = max(dp(i-1,j), dp(i,j-1)) + A(i,j);
            dq(i,j) = min(dq(i-1,j), dq(i,j-1)) + A(i,j);
        end
    end
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
    subplot(1,2,k)
    imagesc(A);
    colormap summer;
    colorbar;
    axis equal tight;
    set(gca,'YDir','normal');
    hold on;
    plot(path_dp(:,2), path_dp(:,1), '-r', 'LineWidth', 2);
    plot(path_dq(:,2), path_dq(:,1), '-b', 'LineWidth', 2);
    legend('Last passage path', 'First passage path', 'Location', 'southoutside');
    
    hold off;
end
sgtitle('Directed First and Last Passage Percolation on Random Grids');
