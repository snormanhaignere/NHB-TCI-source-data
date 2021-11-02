% Data underlying anatomical analyses (Figure 4 and corresponding stats)

% load table with data
% for each electrode, we specify
% the subject ID
% hemisphere
% electrode type: depth (D), surface (S), grid (G)
% distance from PAC in millimeters (only present for electrodes <= 30 mm from PAC)
% the annular bin (equivalent to: floor(distPAC/10)+1)
% the nearest vertex on the FsAverage template brain (which has 163842 total vertices)
% integration window width in seconds
% integration window center in seconds
source_directory = fileparts(fileparts(which('anatomical_data.m')));
T = readtable([source_directory '/data/fig4.xlsx']);

% convert width and center to log scale
T.winwidth = log2(T.winwidth/0.05);
T.wincenter = log2(T.wincenter/0.05);

close all;

%% Linear mixed effects model

stat = {'winwidth', 'wincenter'};
for i = 1:length(stat)
    
    % fit LME model
    % note that electrodes >30 mm have NaN distances and are discarded
    lme_equation = [stat{i} ' ~ electype + distPAC + hemi + (1 + hemi + distPAC | subjid)'];
    lme_model = fitlme(T, lme_equation);
    fprintf('\n\n--------- %s ---------\n\n', stat{i});
    fprintf('%s', 'LME model\n\n%s', evalc('lme_model.disp'));
    
    % F-test for distance
    H = zeros(1, length(lme_model.CoefficientNames));
    H(ismember(lme_model.CoefficientNames, 'distPAC')) = 1;
    [p, F, df1, df2] = coefTest(lme_model, H, 0, 'DFMethod', 'satterthwaite');
    fprintf('\n\n\nF-test for distance: F(%.4f, %.4f)=%.3f, p=%.5f\n', df1, df2, F, p)
    
    % F-test for hemisphere
    H = zeros(1, length(lme_model.CoefficientNames));
    H(ismember(lme_model.CoefficientNames, 'hemi_lh')) = 1;
    [p, F, df1, df2] = coefTest(lme_model, H, 0, 'DFMethod', 'satterthwaite');
    fprintf('F-test for hemisphere: F(%.4f, %.4f)=%.3f, p=%.5f\n', df1, df2, F, p);
    
end

%% Plot annular bin medians with bootstrapped errors, Figure 4b

% statistic to plot (width or center)
stat = T.winwidth;

% unique integer for each subject
[unique_subjids, ~, si] = unique(T.subjid);
[unique_hemis, ~, hi] = unique(T.hemi);

% data matrix
% elec x [subject, hemisphere, annular bin, stat]
D = [si, hi, T.annularbin, stat];

n_bins = 3;
n_smps = 10000;
n_hemi = 2;
bin_smps = nan(n_smps, n_bins, n_hemi);
for j = 1:n_smps
    
    if mod(j,100)==0
        fprintf('Bootstrap %d of %d\n', j, n_smps); drawnow;
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
    
    for i = 1:nbins
        for k = 1:n_hemi
            % all sampled electrodes for a particular bin / subject
            xi = D_bstrap_subj(:,2) == k & D_bstrap_subj(:,3) == i;
            bin_smps(j,i,k) = median(D_bstrap_subj(xi,4));
        end
    end
end

% plot medians and standard error of the sampled distribution
cols = colormap('lines');
cols([1,2],:) = cols([2,1],:);
m = nanmedian(bin_smps,1);
se = stderr_from_samples(bin_smps, 'NaN_frac', 0.1);
figh = figure;
set(figh, 'Position', [200 200 600 450]);
x = (0:2)*10 + 5;
h = nan(1,2);
for q = 1:2
    hold on;
    h(q) = errorbar(x, m(1,:,q), m(1,:,q) - se(1,:,q), se(2,:,q) - m(1,:,q), ...
        '-', 'LineWidth', 2, 'Color', cols(q,:));
    set(h(q), 'Marker', 'none');
end
xlim([0, 30]);
yL = [0.025, 0.4];
ylim(log2(yL/0.05));
yticks = [0.025, 0.05, 0.1, 0.2, 0.4];
set(gca, 'YTick', log2(yticks/0.05), 'YTickLabel', yticks);
legend(h(:), unique_hemis);

%% Center vs. width, Figure 4c

% plot center vs. width (logarithmically transformed)
figure;
plot(T.winwidth, T.wincenter, 'ko'); hold on;

% plot fit, note complications below are due to plotting on log-log scale
cols = lines;
mylog = @(x)(log2(x/0.05));
invlog = @(x)(2.^(x)*0.05);
f = @(x,a,b)(mylog(a*invlog(x) + b));
bounds = mylog([0.025, 0.8]);
x = bounds(1):0.01:bounds(2);
plot(x, f(x, 0.66, 0.021), '-', 'LineWidth', 2, 'Color', cols(2,:)); % best fit
plot(x, f(x, 0.5, 0), '-', 'LineWidth', 2, 'Color', cols(1,:)); % causal minimum
xlim(bounds); ylim(bounds);
set(gca, 'XTick', mylog([50, 100, 200, 400]/1000), 'XTickLabel', [50, 100, 200, 400]);
set(gca, 'YTick', mylog([50, 100, 200, 400]/1000), 'YTickLabel', [50, 100, 200, 400]);
xlabel('Width (ms)'); ylabel('Center (ms)');