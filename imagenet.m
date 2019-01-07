% unittest for Energy distributions

% Version history:
% - 6/8/18: imagenet added to folder 'ImportanceSampling'
% - 1/7/19: clean up file structures & 'imagenet.m' for sharing

str = computer;
% either mac or linux
% use 'mex_this.m' to compile the cpp code for matlab if neither works
if strcmp(str, 'MACI64')
    addpath('source_3D_matlab_mac');
else
    addpath('source_3D_matlab_linux');
end
addpath('source_SMC');
addpath('WARP');

rng(0);

for ith_obs = 1:100
    load(sprintf('imagenet_matlab/obs%d.mat', ith_obs));
    
    obs_true_raw = A;
    obs_true = obs_true_raw;
    
    sigma = 0.1;
    obs = obs_true + randn(size(obs_true, 1)) .* sigma;
    
    
    dir = sprintf('energyDistribution/obs%d', ith_obs);
    mkdir(dir);
    % % hyper0 = hyper_default(obs);
    dimension = size(obs);
    hyper0_t = hyper_default(obs(:), dimension', false);
    n = numel(obs);
    
    
    [lambda_map, C, compact_idx_d, flag_first, flag_not_first, ...
        pointer_flag, rank_left_child, rank_right_child, dividible, ...
        par, y, ss, M_d_map] = Intermediate(obs(:), dimension', hyper0_t);
    
    %% global settings
    total_level = log2(numel(obs));
    num_node = 2 * dimension - 1;
    N = prod(num_node);
    sigma_hat = par(total_level + 1, 3);
    %
    % save('temp.mat')
    
    % Importance sampling
    % periodization:this boundary mode makes c the same length as obs
    dwtmode('per');
    
    % number of random trees
    n_tree = 1;
    
    position_all = cell([n_tree, 1]);
    
    tic
    for ith_tree = 1:n_tree
        % from l = 0 to total_level - 1; last total_level gives the permutation
        node_in = cell([total_level + 1, 1]);
        % node_in contains the 'ranks' of nodes included;
        % (l+1)th cell is for lth level
        node_in{1} = 1;
        % l = total_level: trivial (all 1 by 1 nodes are included)
        % - so do not include
        % the cutting directioin of node_in
        direction_in = cell([total_level, 1]);
        
        lambda_mat = zeros(size(rank_left_child));
        lambda_mat((dividible > 0)) = exp(lambda_map);
        
        for l = 0:(total_level - 1)
            % number of nodes at this level = 2^l;
            rank_node = node_in{l + 1};
            divide_prob = lambda_mat(rank_node, :);
            % if 1, then diretion 1; if 0, then direction 2
            flag = (rand([size(divide_prob, 1), 1]) < divide_prob(:, 1));
            which_direction = 2 - flag;
            % row: rank_node; col: which_direction
            direction_in{l + 1} = which_direction;
            idx_selected = (uint64(which_direction) - 1) * N + rank_node;
            
            rank_node_next_level = cat(1, ...
                rank_left_child(idx_selected)', ...
                rank_right_child(idx_selected)'); % rbind - 2 rows
            node_in{l+2} = rank_node_next_level(:);
        end
        
        
        [i, j] = ind2sub(num_node, node_in{total_level + 1});
        position = sub2ind(size(obs), ...
            i - size(obs, 1) + 1, j - size(obs, 2) + 1);
        % obs_vec = obs(position);\
        position_all{ith_tree} = position;
    end
    
    % summarize one tree
    
    for ith_tree = 1:n_tree
        % ith_tree = 1;
        position = position_all{ith_tree};
        
        % for type = 1:2
        type = 2;
        if type == 1
            input = obs;
        else
            input = obs_true;
        end
        % input = obs;
        % input = obs_true;
        
        for flag_DWT = 1:2
            
            c = cell([2, 1]);
            [c{1}, wave_level] = wavedec(input(position), log2(n), 'db1');
            
            if (flag_DWT == 1)
                [c{2}, ~] = wavedec(input(:), log2(n), 'db1');
            else
                [a,h,v,d] = haart2(input);
                
                c{2} = a;
                for j = 1:(floor(log2(n)/2))
                    c{2} = [c{2}; h{j}(:); v{j}(:); d{j}(:)];
                end
            end
            
            ftrue = figure;
            imshow(A);
            
            % use energy distribution + coefficient saving
            f0 = figure;
            y = zeros([numel(c{1}), 2]);
            for i = 1:2
                y(:, i) = -sort(-abs(c{i}));
                y(:, i) = cumsum(y(:, i).^2);
            end
            
            y = y ./ sum(c{1}.^2);
            x = (1:min([100000, numel(obs_true)]));
            energy_grid = linspace(0.85, 0.95, 500);
            vq1 = rb_interp1(y(x, 1),x, energy_grid);
            vq2 = rb_interp1(y(x, 2),x, energy_grid);
            yyaxis right;
            plot(energy_grid,vq1, '--', 'color', 'blue', 'Linewidth', 2);
            hold on
            plot(energy_grid,vq2, '.-', 'color', 'red', 'Linewidth', 2);
            hold off
            ylabel('Number of coefficients')
            
            yyaxis left;
            plot(energy_grid, 1 - vq1 ./ vq2, '-', 'color', 'black', 'LineWidth', 3);
            ylabel('Percentage of coefficients saving');
            a=cellstr(num2str(get(gca,'ytick')'*100));
            pct = char(ones(size(a,1),1)*'%');
            new_yticks = [char(a),pct];
            set(gca,'yticklabel',new_yticks)
            
            %             % use energy distribution
            %
            %             y = zeros([numel(c{1}), 2]);
            %             for i = 1:2
            %                 y(:, i) = -sort(-abs(c{i}));
            %                 y(:, i) = cumsum(y(:, i).^2);
            %             end
            %
            %             y = y ./ sum(c{1}.^2);
            %
            %             f1 = figure;
            %             x = (1:1000);
            %             plot(x,y(x, 1), 'color', 'blue', 'Linewidth', 2);
            %             hold on
            %             plot(x,y(x, 2), 'color', 'red', 'Linewidth', 2)
            %             hold off
            %             % legend('WARP', 'Determinstic Vectorization')
            %
            %             % plot the raw coefficients
            %             y = zeros([numel(c{1})-1, 2]);
            %             for i = 1:2
            %                 y(:, i) = sort(c{i}(2:end),'descend','ComparisonMethod','abs');
            %             end
            %             f2 = figure;
            %             x = (2:1000);
            %             % x = (2:(numel(c{1}) - 1));
            %             plot(x,y(x, 2), '.', 'Color', 'red'); yl = ylim;
            %             % plot(x,y(x, 2), '.', 'Color', [137 137 137]./256); yl = ylim;
            %             hold on;
            %             plot(x,y(x, 1), '.', 'Color', 'blue');
            %             hold off;
            
            % use sparsity
            
            threshold = sigma * sqrt(2 * log(n));
            
            % x_pt = [0.01, 0.05, 0.1, 0.5, 1];
            x_pt = (10:100)./200;
            y = zeros([numel(x_pt), 2]);
            for k = 1:numel(x_pt)
                for i = 1:2
                    y(k,i) = mean(abs(c{i}) < threshold * x_pt(k));
                end
            end
            
            %             f3 = figure;
            %             plot(x_pt,y(:, 1),'color', 'blue', 'Linewidth', 2);
            %             hold on
            %             plot(x_pt,y(:, 2),'color', 'red', 'Linewidth', 2)
            %             % legend('WARP', 'Determinstic Vectorization')
            %             hold off;
            saveas(ftrue, sprintf('%s/true%dtype%d%dd', ...
                dir,ith_tree, type, flag_DWT), 'png');
            
            saveas(f0, sprintf('%s/saving%dtype%d%dd', ...
                dir,ith_tree, type, flag_DWT), 'png');
            %             saveas(f1, sprintf('%s/energy%dtype%d%dd', ...
            %                 dir,ith_tree, type, flag_DWT), 'png');
            %             saveas(f2, sprintf('%s/coef%dtype%d%dd', ...
            %                 dir, ith_tree, type, flag_DWT), 'png');
            %             saveas(f3, sprintf('%s/sparsity%dtype%d%dd', ...
            %                 dir, ith_tree, type, flag_DWT), 'png');
            
            close all
        end
    end
end

