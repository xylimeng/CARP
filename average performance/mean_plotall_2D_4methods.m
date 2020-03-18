
clear all
load('cr_2.mat')
load('cr_4.mat')
load('cr_5.mat')
load('cr_6.mat')

load('psnr_2.mat')
load('psnr_4.mat')
load('psnr_5.mat')
load('psnr_6.mat')
%2.e2e;4.CARP;5.JPEG;6.JPEG2000
%% find all index need to keep
% method 1
ind_1 = find(min(cr_4, [], 2)<15);
ind_2 = find(max(cr_4, [], 2)>45);
ind_tmp0 = intersect(ind_1,ind_2);

ind_1 = find(min(cr_5, [], 2)<15);
ind_2 = find(max(cr_5, [], 2)>45);
ind_tmp1 = intersect(ind_1,ind_2);

ind = intersect(ind_tmp0, ind_tmp1);

% method 2
% ind_1 = find(min(cr_2, [], 2)<10);
% ind_2 = find(max(cr_2, [], 2)>50);
% ind = intersect(ind_1,ind_2);
cr_2_new = cr_2(ind, :);
psnr_2_new = psnr_2(ind, :);
cr_all = cr_2_new(:);
cr_min = max(min(cr_2_new, [], 2));
cr_max = min(max(cr_2_new, [], 2));
cr_2_share = sort(cr_all(intersect(find(cr_all>=cr_min), find(cr_all<=cr_max))));
psnr_2_share = zeros(length(ind), length(cr_2_share));
for k=1:length(ind)
    [cr_2_s,idx]=sort(cr_2_new(k,:));
    psnr_2_s = psnr_2_new(k,idx);
    psnr_2_share(k,:) = linear_fit(cr_2_share, cr_2_s, psnr_2_s);
end
psnr_2_mean = mean(psnr_2_share);

% method 4
% ind_1 = find(min(cr_4, [], 2)<10);
% ind_2 = find(max(cr_4, [], 2)>50);
% ind = intersect(ind_1,ind_2);
cr_4_new = cr_4(ind, :);
psnr_4_new = psnr_4(ind, :);
cr_all = cr_4_new(:);
cr_min = max(min(cr_4_new, [], 2));
cr_max = min(max(cr_4_new, [], 2));
cr_4_share = sort(cr_all(intersect(find(cr_all>=cr_min), find(cr_all<=cr_max))));
psnr_4_share = zeros(length(ind), length(cr_4_share));
for k=1:length(ind)
    [cr_4_s,idx]=sort(cr_4_new(k,:));
    psnr_4_s = psnr_4_new(k,idx);
    psnr_4_share(k,:) = linear_fit(cr_4_share, cr_4_s, psnr_4_s);
end
psnr_4_mean = mean(psnr_4_share);
% method 5
% ind_1 = find(min(cr_5, [], 2)<10);
% ind_2 = find(max(cr_5, [], 2)>50);
% ind = intersect(ind_1,ind_2);
cr_5_new = cr_5(ind, :);
psnr_5_new = psnr_5(ind, :);
cr_all = cr_5_new(:);
cr_min = max(min(cr_5_new, [], 2));
cr_max = min(max(cr_5_new, [], 2));
cr_5_share = sort(cr_all(intersect(find(cr_all>=cr_min), find(cr_all<=cr_max))));
psnr_5_share = zeros(length(ind), length(cr_5_share));
for k=1:length(ind)
    [cr_5_s,idx]=sort(cr_5_new(k,:));
    psnr_5_s = psnr_5_new(k,idx);
    psnr_5_share(k,:) = linear_fit(cr_5_share, cr_5_s, psnr_5_s);
end
psnr_5_mean = mean(psnr_5_share);
% method 6
% ind_1 = find(min(cr_6, [], 2)<10);
% ind_2 = find(max(cr_6, [], 2)>50);
% ind = intersect(ind_1,ind_2);
cr_6_new = cr_6(ind, :);
psnr_6_new = psnr_6(ind, :);
cr_all = cr_6_new(:);
cr_min = max(min(cr_6_new, [], 2));
cr_max = min(max(cr_6_new, [], 2));
cr_6_share = sort(cr_all(intersect(find(cr_all>=cr_min), find(cr_all<=cr_max))));
psnr_6_share = zeros(length(ind), length(cr_6_share));
for k=1:length(ind)
    [cr_6_s,idx]=sort(cr_6_new(k,:));
    psnr_6_s = psnr_6_new(k,idx);
    psnr_6_share(k,:) = linear_fit(cr_6_share, cr_6_s, psnr_6_s);
end
psnr_6_mean = mean(psnr_6_share);

% plot
plot(cr_4_share, psnr_4_mean, cr_5_share, psnr_5_mean, ...
    cr_6_share, psnr_6_mean, cr_2_share, psnr_2_mean)
xlim([0 50]) 
xlabel('Compression Ratio'); ylabel('PSNR');
legend('CARP','JPEG','JPEG2000','End-to-End')
