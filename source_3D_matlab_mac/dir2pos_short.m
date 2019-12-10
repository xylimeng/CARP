function [position, prun_block]=dir2pos_short(direction, prun_num, dimension)
prun_block = cell(prun_num, 1);
m = length(dimension);
J = log2(prod(dimension));
level = J+1;
dim_level_old = cell(1,1);
dim_level_old{1,1} = dimension;
ind = reshape(1:prod(dimension), dimension);
ind_level_old = cell(1,1);
ind_level_old{1,1} = ind;
dir_l = direction(1);
node_in = 1;
long_prun = 0;
prun_k = 1;
for l=1:(level-1)
    fprintf('Total level: %i. Current level: %i.\n',level,l);
    %dir_l = direction(2^(l-1):(2^l-1));
    len_level = length(dir_l);
    long_dir_new = -2*ones(1,2*len_level);
    long_prun_new = zeros(1,2*len_level);
    dim_level_new = cell(1,2*len_level);
    ind_level_new = cell(1,2*len_level);
    for t=1:len_level
        dim0 = dim_level_old{1,t};
        ind0 = ind_level_old{1,t};
        dir0 = dir_l(t)+1;
        if dir0 >0
            cut_point = floor(dim0(dir0)/2);
            dim_tmp = dim0;
            dim_tmp(dir0) = cut_point;
            dim_level_new{1,2*(t-1)+1} = dim_tmp;
            dim_level_new{1,2*t} = dim_tmp;
        end
        if m==2
            if dir0==1
                ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :);
                ind_level_new{1,2*t}= ind0((cut_point+1):end, :);
            elseif dir0==2
                ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point);
                ind_level_new{1,2*t}= ind0(:, (cut_point+1):end);
            elseif dir0==0
                if dim0(2)>1
                    cut_point = floor(dim0(2)/2);
                    dim_tmp = dim0;
                    dim_tmp(2) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point);
                    ind_level_new{1,2*t} = ind0(:, (cut_point+1):end);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(1)>1
                    cut_point = floor(dim0(1)/2);
                    dim_tmp = dim0;
                    dim_tmp(1) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :);
                    ind_level_new{1,2*t} = ind0((cut_point+1):end, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                end
                if l>1 && long_prun(t)==1
                    prun_block{prun_k, 1} = ind0;
                    prun_k = prun_k+1;
                end
            end
        elseif m==3
            if dir0==1
                ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :, :);
                ind_level_new{1,2*t}= ind0((cut_point+1):end, :, :);
            elseif dir0==2
                ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point, :);
                ind_level_new{1,2*t}= ind0(:, (cut_point+1):end, :);
            elseif dir0==3
                ind_level_new{1,2*(t-1)+1}= ind0(:, :, 1:cut_point);
                ind_level_new{1,2*t}= ind0(:, :, (cut_point+1):end);
            elseif dir0==0
                if dim0(3)>1
                    cut_point = floor(dim0(3)/2);
                    dim_tmp = dim0;
                    dim_tmp(3) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, :, 1:cut_point);
                    ind_level_new{1,2*t} = ind0(:, :, (cut_point+1):end);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(2)>1
                    cut_point = floor(dim0(2)/2);
                    dim_tmp = dim0;
                    dim_tmp(2) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point, :);
                    ind_level_new{1,2*t} = ind0(:, (cut_point+1):end, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(1)>1
                    cut_point = floor(dim0(1)/2);
                    dim_tmp = dim0;
                    dim_tmp(1) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :, :);
                    ind_level_new{1,2*t} = ind0((cut_point+1):end, :, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                end
                if l>1 && long_prun(t)==1
                    prun_block{prun_k, 1} = ind0;
                    prun_k = prun_k+1;
                end
            end
        elseif m==4
            if dir0==1
                ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :, :, :);
                ind_level_new{1,2*t}= ind0((cut_point+1):end, :, :, :);
            elseif dir0==2
                ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point, :, :);
                ind_level_new{1,2*t}= ind0(:, (cut_point+1):end, :, :);
            elseif dir0==3
                ind_level_new{1,2*(t-1)+1}= ind0(:, :, 1:cut_point, :);
                ind_level_new{1,2*t}= ind0(:, :, (cut_point+1):end, :);
            elseif dir0==4
                ind_level_new{1,2*(t-1)+1}= ind0(:, :, :, 1:cut_point);
                ind_level_new{1,2*t}= ind0(:, :, :, (cut_point+1):end);
            elseif dir0==0
                if dim0(4)>1
                    cut_point = floor(dim0(4)/2);
                    dim_tmp = dim0;
                    dim_tmp(4) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, :, :, 1:cut_point);
                    ind_level_new{1,2*t} = ind0(:, :, :, (cut_point+1):end);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(3)>1
                    cut_point = floor(dim0(3)/2);
                    dim_tmp = dim0;
                    dim_tmp(3) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, :, 1:cut_point, :);
                    ind_level_new{1,2*t} = ind0(:, :, (cut_point+1):end, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(2)>1
                    cut_point = floor(dim0(2)/2);
                    dim_tmp = dim0;
                    dim_tmp(2) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(:, 1:cut_point, :, :);
                    ind_level_new{1,2*t} = ind0(:, (cut_point+1):end, :, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                elseif dim0(1)>1
                    cut_point = floor(dim0(1)/2);
                    dim_tmp = dim0;
                    dim_tmp(1) = cut_point;
                    dim_level_new{1,2*(t-1)+1} = dim_tmp;
                    dim_level_new{1,2*t} = dim_tmp;
                    ind_level_new{1,2*(t-1)+1}= ind0(1:cut_point, :, :, :);
                    ind_level_new{1,2*t} = ind0((cut_point+1):end, :, :, :);
                    long_dir_new((2*(t-1)+1):(2*t)) = -1;
                end
                if l>1 && long_prun(t)==1
                    prun_block{prun_k, 1} = ind0;
                    prun_k = prun_k+1;
                end
            end
        end
    end
    if l<(level-1)
        empty_indx = find(long_dir_new==-2);
        %length(empty_indx)
        if isempty(empty_indx)==0
            empty_len = length(empty_indx);
            long_dir_new(empty_indx)=direction((node_in+1):(node_in+empty_len));
            long_prun_new(empty_indx)=1*(direction((node_in+1):(node_in+empty_len))==-1);
            dir_l = long_dir_new;
            long_prun = long_prun_new;
            node_in = node_in+empty_len;
        else
            dir_l = long_dir_new;
            long_prun = long_prun_new;
            node_in = node_in+0;
        end
    end
    dim_level_old = dim_level_new;
    ind_level_old = ind_level_new;
end
position = cell2mat(ind_level_new);

