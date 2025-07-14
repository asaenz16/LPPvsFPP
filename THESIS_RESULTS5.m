% Conditional Data PLOTS (end‐point and third‐point)
clear;
n    = 300;
m    = n/3;       % third‐point index (100)
t    = 100;

% Preallocate storage
dp_end   = zeros(t,1);
dq_end   = zeros(t,1);
dp_mid   = zeros(t,1);
dq_mid   = zeros(t,1);

% 1) Run trials and record both final and mid‐point passage times
for k = 1:t
    A = exprnd(1, [n,n]);
    
    % MAX path
    dp = zeros(n);
    dp(1,:) = cumsum(A(1,:));
    dp(:,1) = cumsum(A(:,1));
    for i = 2:n
        for j = 2:n
            dp(i,j) = max(dp(i,j-1), dp(i-1,j)) + A(i,j);
        end
    end
    
    % MIN path
    dq = zeros(n);
    dq(1,:) = cumsum(A(1,:));
    dq(:,1) = cumsum(A(:,1));
    for i = 2:n
        for j = 2:n
            dq(i,j) = min(dq(i,j-1), dq(i-1,j)) + A(i,j);
        end
    end
    
    dp_end(k) = dp(n,n);
    dq_end(k) = dq(n,n);
    dp_mid(k) = dp(m,m);      % record at (150,150)
    dq_mid(k) = dq(m,m);
end

% 2) Discretize into integer bins
idx_Le = round(dp_end);
idx_Fe = round(dq_end);
idx_Lm = round(dp_mid);
idx_Fm = round(dq_mid);

% 3) Joint histograms → pmfs
maxLe = max(idx_Le); maxFe = max(idx_Fe);
maxLm = max(idx_Lm); maxFm = max(idx_Fm);

pmf_e = accumarray([idx_Fe, idx_Le], 1, [maxFe, maxLe]) / t;
pmf_m = accumarray([idx_Fm, idx_Lm], 1, [maxFm, maxLm]) / t;

% 4) Conditional P(L|F) for each
row_e = sum(pmf_e,2);
cond_e = pmf_e ./ row_e; cond_e(row_e==0,:) = 0;

row_m = sum(pmf_m,2);
cond_m = pmf_m ./ row_m; cond_m(row_m==0,:) = 0;

% 5) Side‐by‐side heatmaps
figure;

% -- N = 300 endpoint --
subplot(1,2,2);
imagesc(1:maxFe, 1:maxLe, cond_e');
set(gca,'YDir','normal');
set(gca,'CLim',[0 0.15]);
xlabel('First Passage Time');
ylabel('Last Passage Time');
title('Conditional PMF (N=300)');
colorbar;
axis([min(idx_Fe) max(idx_Fe) min(idx_Le) max(idx_Le)]);

% -- N = 100  --
subplot(1,2,1);
imagesc(1:maxFm, 1:maxLm, cond_m');
set(gca,'YDir','normal');
set(gca,'CLim',[0 0.2]);
xlabel('First Passage Time');
ylabel('Last Passage Time');
title('Conditional PMF (N=100)');
colorbar;
axis([min(idx_Fm) max(idx_Fm) min(idx_Lm) max(idx_Lm)]);
