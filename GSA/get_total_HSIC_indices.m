function [T_hsic] = get_total_HSIC_indices(generate_X, evaluate_Y, N_samples, N_inputs, bandwidth)
% GET_TOTAL_HSIC_INDICES Computes Total HSIC indices with O(N^2) memory footprint.
%
% INPUTS:
%   generate_X - Function handle @(N) returning N x d inputs.
%   evaluate_Y - Function handle @(X) returning N x k outputs.
%   N_samples  - Integer, number of Monte Carlo samples.

    % --- 1. Generate Data ---

%    fprintf('Generating %d input samples...\n', N_samples);
    X = generate_X(N_samples, N_inputs);
    [m, n] = size(X);
    
%    fprintf('Evaluating model outputs...\n');
    Y = evaluate_Y(X);
    
    J = (1/m) * ones(m,m);
    
    % --- 2. Precompute the Output Kernel (L_final) ---
%    fprintf('Computing Output Kernel...\n');
    DY = pdist2(Y, Y);
    upperY = DY(triu(true(m),1));
    
    if bandwidth ==1
        sigY = var(Y);
    elseif bandwidth ==0
        sigY = max(median(upperY(:)), 1e-8);
    end
    
    L_y_gaussian = exp(-(DY.^2) / (2 * sigY^2));
    clear DY upperY; % Free RAM immediately
    
    L_y_centered = L_y_gaussian - J*L_y_gaussian - L_y_gaussian*J + J*L_y_gaussian*J;
    clear L_y_gaussian;
    
    L_final = L_y_centered + ones(m,m);
    clear L_y_centered;
    
    % --- 3. Precompute Input Bandwidths (sigX) ---
%    fprintf('Precomputing input bandwidths...\n');
    sigX = zeros(1, n); 
    for d = 1:n
        if bandwidth == 1
            sigX(d) = sqrt(n*var(X(:,d)));

        elseif bandwidth ==0
            Dx_d = pdist2(X(:,d), X(:,d));
            upper_d = Dx_d(triu(true(size(Dx_d)),1));
            sigX(d) = median(upper_d(:));
        end
    end
    clear Dx_d upper_d; 
    
    % --- 4. Build Total Input Gram Matrix (Pass 1) ---
%    fprintf('Building Total Input Gram Matrix on the fly...\n');
    K = ones(m,m); 
    for d = 1:n
        % Compute single dimension kernel
        Dx_d = pdist2(X(:,d), X(:,d));
        K_gauss = exp(-(Dx_d.^2) / (2 * sigX(d)^2)); %Gaussian
        %K_gauss = exp(-(Dx_d).^(1/2) / (2 * sigX(d)^2)); %Laplace
        
        % Center and add 1
        K_cent = K_gauss - J*K_gauss - K_gauss*J + J*K_gauss*J;
        K_final_d = K_cent + ones(m,m);
        
        % Multiply into the running total
        K = K .* K_final_d; 
    end
    clear Dx_d K_gauss K_cent K_final_d;
    
    % --- 5. Compute Total HSIC ---
%    fprintf('Calculating Total and Partial HSIC indices...\n');
    HSIC_XY = compute_HSIC_internal(K, L_final);
    
    % --- 6. Compute Partial HSICs (Pass 2) ---
    T_hsic = zeros(1, n); 
    parfor i = 1:n 
        % Re-compute just the i-th kernel on the fly
        Dx_i = pdist2(X(:,i), X(:,i));
        K_gauss = exp(-(Dx_i.^2) / (2 * sigX(i)^2)); %Gaussian
        %K_gauss = exp(-(Dx_i).^(1/2) / (2 * sigX(i)^2)); %Laplacian
        K_cent = K_gauss - J*K_gauss - K_gauss*J + J*K_gauss*J;
        K_i = K_cent + ones(m,m);
        
        % Divide out the i-th component
        K_noti = K ./ (K_i+1e-12); 
        
        % Compute HSIC for this subset
        HSIC_Xnoti = compute_HSIC_internal(K_noti, L_final);
        T_hsic(i) = 1 - HSIC_Xnoti / HSIC_XY;
    end
%    fprintf('HSIC computation complete.\n');
end

% --- Internal Helper Function ---
function [HSIC] = compute_HSIC_internal(K, L)
    [n,~] = size(K);
    Lsum = sum(L, 1); 
    
    HSIC = 0;
    for i=1:n
        HSIC = HSIC + dot(K(i,:)', L(:,i)) - (2/n)*sum(K(:,i))*sum(L(:,i)) + (1/n^2)*sum(K(:,i))*sum(Lsum);
    end
    HSIC = HSIC / n^2;
end