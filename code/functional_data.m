% Data underlying functional analyses (Figure 5 and corresponding stats)

close all;

%% PC scatter, Figure 5b

% load and format data
T = readtable([source_directory '/data/fig5_PCs.xlsx']);
n_groups = 3;
n_PCs = 2;
n_sounds = 119; 
pc_data = [T.group1_pc1, T.group1_pc2, T.group2_pc1, T.group2_pc2, T.group3_pc1, T.group3_pc2];
pc_data = reshape(pc_data, [n_sounds, n_PCs, n_groups]);

% plot
category_labels = {'drumming', 'music', 'song', 'englishspeech', ...
    'foreignspeech', 'nonspeechvoc', 'anivoc', 'aninonvoc', 'nature', ...
    'mechanical', 'envsounds'};
colors = [...
         0    0.7255    0.8980; ...
    0.0784    0.1686    0.5490; ...
    1.0000    0.2000    0.2000; ...
    0.0549    0.3294    0.2353; ...
    0.4039    0.7059    0.5686; ...
    0.3529    0.1804    0.5137; ...
    0.5961    0.3765    0.8157; ...
    0.9059    0.3451    0.4510; ...
    0.9216    0.9216    0.2863; ...
    0.8706    0.4902         0; ...
    0.3922    0.3922    0.3922; ...
    ];
figh = figure;
set(figh, 'Position', [100, 100, 800, 250]);
group_names = {'short', 'intermediate', 'long'};
for i = 1:n_groups
    subplot(1, 3, i);
    for k = 1:n_sounds
        plot(pc_data(k,1,i), pc_data(k,2,i), 'o', ...
            'Color', colors(T.category_index(k),:), 'LineWidth', 2); hold on;
    end
    xlabel('PC1'); ylabel('PC2');
    title(group_names{i});
end

%% Linear mixed effects model for prediction accuracies

% load prediction data
source_directory = fileparts(fileparts(which('anatomical_data.m')));
T = readtable([source_directory '/data/fig5_predictions.xlsx']);

% convert width and center to log scale
T.winwidth = log2(T.winwidth/0.05);
T.wincenter = log2(T.wincenter/0.05);

% estimate LME model
lme_model = fitlme(T, 'rcatminuscoch ~ winwidth + electype + hemi + (1 + winwidth + hemi | subjid)');

% determine significance for distance to PAC
H = zeros(1, length(lme_model.CoefficientNames));
H(ismember(lme_model.CoefficientNames, 'winwidth')) = 1;
[p, F, df1, df2] = coefTest(lme_model, H, 0, 'DFMethod', 'satterthwaite');
fprintf('\n\n\nF-test for distance: F(%.4f, %.4f)=%.3f, p=%.5f\n\n\n', df1, df2, F, p);

%% Bootstrapped, noise-corrected r^2 for each electrode group, Figure 5c

sign_and_square = @(x)sign(x).*x.^2;
stats = sign_and_square([T.rcoch, T.rcat, T.rcochcat]./T.rceil);

% unique integer for each subject
[unique_subjids, ~, si] = unique(T.subjid);

% data matrix
% elec x [subject, group, statistics]
D = [si, T.group, stats];

% bootstrapped medians for each groups
n_groups = max(T.group);
n_stats = size(stats,2);
n_smps = 10000;
group_medians_bstrap_smps = nan(n_smps, n_groups, n_stats);
for j = 1:n_smps
    
    if mod(j,100)==0
        % fprintf('Bootstrap %d of %d\n', j, n_smps); drawnow;
    end
    
    % bootstrap electrodes separately for each subject each subject
    D_bstrap_elec = [];
    for i = 1:max(si)
        % resample all electrodes from a single subject
        elecs_from_single_subj = find(D(:,1) == i);
        resampled_elecs = randsample(elecs_from_single_subj, length(elecs_from_single_subj), 1);
        
        % save the statistic, bin, and subject index for these electrodes
        D_bstrap_elec = cat(1, D_bstrap_elec, D(resampled_elecs, :));
    end
    
    % bootstrap subjects, replicating the resampled electrodes
    D_bstrap_subj = [];
    resampled_subjects = randsample(max(si), max(si), 1); % sample subjects with replacement
    for i = 1:max(si)
        % all of the previously resample electrodes from a single subject
        elecs_from_single_subj = D_bstrap_elec(:,1) == resampled_subjects(i);
        
        % save the resampled electrodes from the sampled subject
        D_bstrap_subj = cat(1, D_bstrap_subj, D_bstrap_elec(elecs_from_single_subj, :));
    end
    
    for k = 1:n_groups
        % all sampled electrodes for a particular group
        xi = D_bstrap_subj(:,2) == k;
        group_medians_bstrap_smps(j,k,:) = median(D_bstrap_subj(xi,3:end));
    end
end

% plot medians and standard errors across samples
m = nanmedian(group_medians_bstrap_smps);
se = stderr_from_samples(group_medians_bstrap_smps, 'NaN_frac', 0.1);
cmap = colormap('lines');
figh = figure;
set(figh, 'Position', [100 100 300 300]);
hold on;
h = nan(n_groups, n_stats);
for j = 1:n_groups
    for k = 1:n_stats
        h(j,k) = bar(k + (j-1)*(n_stats+1), m(1,j,k), 'FaceColor', cmap(k,:));
        he = errorbar(k + (j-1)*(n_stats+1), m(1,j,k), m(1,j,k)-se(1,j,k), se(2,j,k)-m(1,j,k), 'k.', 'LineWidth', 2);
        set(he, 'Marker', 'none');
    end
end
xlim([0, n_groups*(n_stats+1)]); ylim([0, 1]);
set(gca, 'XTick', [2, 6, 10], 'XTickLabel', {'Short', 'Intermed', 'Long'});

% print medians and confidence intervals
CI = central_interval_from_samples(group_medians_bstrap_smps, 0.05, 'NaN_frac', 0.1);
group_names = {'short', 'intermediate', 'long'};
stat_names = {'coch', 'cat', 'coch+cat'};
for j = 1:n_groups
    for i = 1:n_stats
        fprintf('%s, %s: median=%.2f, CI=[%.3f, %.3f]\n', ...
            group_names{j}, stat_names{i}, m(1, j, i), CI(1,j,i), CI(2,j,i));
    end
end

%% Scatter plot of unique variance, Fig 5d

% unique variance explained by cochleagram features and category labels
sign_and_square = @(x)sign(x).*x.^2;
unique_coch = sign_and_square(T.rcochcat./T.rceil) - sign_and_square(T.rcat./T.rceil);
unique_cat = sign_and_square(T.rcochcat./T.rceil) - sign_and_square(T.rcoch./T.rceil);

% plot scatter
figure;
h = nan(1,2);
cols = colormap('lines');
h(1) = plot(T.winwidth, unique_coch, 'o', 'Color', cols(1,:), 'LineWidth', 2);
hold on;
h(2) = plot(T.winwidth, unique_cat, 'o', 'Color', cols(2,:), 'LineWidth', 2);

% plot logistic fits
flogist = @(x,a,b,c)(c./(1+exp(-b.*(x-a))));
fcoch = @(x)flogist(x, 1.998, -4.601, 0.206);
fcat = @(x)flogist(x, 2.011, 4.125, 0.332);
x = min(T.winwidth):0.01:max(T.winwidth);
plot(x, fcoch(x), '-', 'Color', cols(1,:), 'LineWidth', 2);
plot(x, fcat(x), '-', 'Color', cols(2,:), 'LineWidth', 2);
set(gca, 'XTick', [0, 1, 2, 3], 'XTickLabel', 50.*2.^([0, 1, 2, 3]));
ylabel('Unique variance');
xlabel('Integration width');
legend(h, {'Coch', 'Cat'});