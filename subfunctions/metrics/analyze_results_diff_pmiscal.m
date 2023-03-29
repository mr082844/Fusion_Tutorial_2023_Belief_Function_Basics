function analyze_results_diff_pmiscal(varargin)
% read in variables
optin = {'./tests/results_rand_geo_perfect_pmiscal.mat'};
if 0 < nargin
    optin(1:nargin) = varargin(:);
elseif 1 < nargin
    return;
end
[load_this_file] = optin{:};

load(load_this_file);
field_names = fields(results);
idx_remove = contains(field_names,'init');
field_names(idx_remove) = [];
N_src_list = unique(results.(field_names{1}).N_src);
N_src_list = N_src_list(~isnan(N_src_list));
delta_p = 0.01;
prob_bins = 0 : delta_p : 1;
Pmc_list = flip(unique(results.(field_names{1}).Pmiscal));
Pmc_list = Pmc_list(~isnan(Pmc_list));
for ialg = 1 : length(field_names)
    alg_name = field_names{ialg};
    zeros_init = zeros(length(N_src_list),1);
    ones_init = ones(length(N_src_list),1);
    cell_init = cell(length(N_src_list),1);
    DKL.(alg_name).DS = zeros_init;
    DKL.(alg_name).REL = zeros_init;
    DKL.(alg_name).RES = zeros_init;
    DKL.(alg_name).UNC = ones_init;
    DKL.(alg_name).confMat = cell_init;
    DKL.(alg_name).error = zeros_init;
    DKL.(alg_name).noclass = zeros_init;
    DKL.(alg_name).correct = zeros_init;
    DKL.(alg_name).PCAL = zeros_init;
    DKL.(alg_name).PCAL2 = zeros_init;
    DKL.(alg_name).avg_delta = zeros_init;
    DKL.(alg_name).var_delta = zeros_init;
    DKL.(alg_name).obs = cell_init;
    DKL.(alg_name).obswt = cell_init;
    
    for iNS = 1:length(N_src_list)
        N_sources = N_src_list(iNS);
        idx_rows = results.(alg_name).N_src == N_sources ...
            & ~isnan(results.(alg_name).fused_prob(:,1))...
            & ~isnan(results.(alg_name).fused_prob(:,2));
        N_rows.(alg_name) = sum(idx_rows);
        
        prob = results.(alg_name).fused_prob(idx_rows,1);
        pred = zeros(size(prob));
        pred(prob>0.51) = 1;
        pred(prob<0.49) = 2;
        DKL.(alg_name).confMat{iNS} = ...
            confusionmat(results.(alg_name).truth(idx_rows)...
            ,pred);
        if ~isempty(DKL.(alg_name).confMat{iNS})
            DKL.(alg_name).error(iNS) = (DKL.(alg_name).confMat{iNS}(end-1,end) ...
                + DKL.(alg_name).confMat{iNS}(end,end-1)) ...
                / sum(sum(DKL.(alg_name).confMat{iNS}));
            DKL.(alg_name).correct(iNS) = (DKL.(alg_name).confMat{iNS}(end-1,end-1) ...
                + DKL.(alg_name).confMat{iNS}(end,end)) ...
                / sum(sum(DKL.(alg_name).confMat{iNS}));
            DKL.(alg_name).noclass(iNS) = 1 ...
                - DKL.(alg_name).correct(iNS) ...
                - DKL.(alg_name).error(iNS);
        end
        for ir = find(idx_rows')
            truth_temp = false(1,2);
            truth_temp(results.(alg_name).truth(ir)) = true;
            DKL.(alg_name).DS(iNS) = DKL.(alg_name).DS(iNS) ...
                + KLD(results.(alg_name).fused_prob(ir,:)...
                ,truth_temp);
        end
        DKL.(alg_name).DS(iNS) = DKL.(alg_name).DS(iNS) ...
            / N_rows.(alg_name);
        
        MC_pc_list = unique(results.(alg_name).MC_pc(idx_rows,1));
        MC_t_list = unique(results.(alg_name).MC_t(idx_rows,1));
        delta = NaN(length(MC_pc_list),length(MC_t_list),2);
        for iMC_pc = 1:length(MC_pc_list)
            MC_pc = MC_pc_list(iMC_pc);
            idx_MC_pc = results.(alg_name).MC_pc == MC_pc;
            for iMC_t = 1 : length(MC_t_list)
                MC_t = MC_t_list(iMC_t);
                idx_MC_t = results.(alg_name).MC_t == MC_t;
                for ic = 1 : 2
                    idx_c = results.(alg_name).truth == ic;
                    idx_test = idx_MC_pc & idx_MC_t & idx_rows & idx_c;
                    if any(idx_test)
                        delta(iMC_pc,iMC_t) = max(results.(alg_name).fused_prob(idx_test,1)) ...
                            - min(results.(alg_name).fused_prob(idx_test,1));
                    end
                end
            end
        end
        delta = reshape(delta,1,numel(delta));
        DKL.(alg_name).avg_delta(iNS) = mean(delta,'omitnan');
        DKL.(alg_name).var_delta(iNS) = var(delta,'omitnan');
        
        DKL.(alg_name).ok_counts = zeros(2,length(0 : delta_p : 1));
        
        idx_rows = results.(alg_name).N_src == N_sources;
        N_rows.(alg_name) = sum(idx_rows);
        
        for ir = find(idx_rows')
            for k1 = 1 : length(prob_bins)
                k2 = length(prob_bins) - k1 + 1;
                pk1 = prob_bins(k1);
                if results.(alg_name).fused_prob(ir,1) ...
                        > pk1 - delta_p/2 ...
                        && results.(alg_name).fused_prob(ir,1) ...
                        < pk1 + delta_p/2
                    DKL.(alg_name).ok_counts(results.(alg_name).truth(ir),k1) = ...
                        DKL.(alg_name).ok_counts(results.(alg_name).truth(ir),k1) + 1;
                    break;
                end
            end
        end
        DKL.(alg_name).ok = zeros(size(DKL.(alg_name).ok_counts));
        DKL.(alg_name).ok = DKL.(alg_name).ok_counts ...
            ./ repmat(sum(DKL.(alg_name).ok_counts,1),2,1);
        DKL.(alg_name).ok_weight = sum(DKL.(alg_name).ok_counts,1) ...
            / sum(sum(DKL.(alg_name).ok_counts));
       
        
        obs1 = zeros(size(DKL.(alg_name).ok(1,:)));
        obs_wgt1 = zeros(size(DKL.(alg_name).ok(1,:)));
        for io = 1 : length(obs1)
            iomn = max(1,io-2);
            iomx = min(length(obs1),io+2);
            obs_wgt1(io) = sum(DKL.(alg_name).ok_weight(iomn:iomx));
            obs1(io) = sum( DKL.(alg_name).ok(1,iomn:iomx)...
                .* DKL.(alg_name).ok_weight(iomn:iomx) ) ...
                / sum(DKL.(alg_name).ok_weight(iomn:iomx));
        end
        temp_wt = flip(DKL.(alg_name).ok_weight);
        obs2 = zeros(size(DKL.(alg_name).ok(1,:)));
        obs_wgt2 = zeros(size(DKL.(alg_name).ok(1,:)));
        for io = length(obs2) : -1 : 1
            iomn = max(1,io-2);
            iomx = min(length(obs1),io+2);
            obs_wgt2(io) = sum(temp_wt(iomn:iomx));
            obs2(io) = sum( DKL.(alg_name).ok(1,iomn:iomx)...
                .* temp_wt(iomn:iomx) ) ...
                / sum(temp_wt(iomn:iomx));
        end
        
        obs = (obs1.*obs_wgt1 + obs2.*obs_wgt2) ./ sum([obs_wgt2;obs_wgt1],1);
        obs_wgt = (obs_wgt1 + obs_wgt2);
        obs2 = zeros(size(DKL.(alg_name).ok(1,:)));
        obs_wgt2 = zeros(size(DKL.(alg_name).ok(1,:)));
        for io = 1 : length(obs)
            iomn = max(1,io-2);
            iomx = min(length(obs),io+2);
            obs_wgt2(io) = sum(obs_wgt(iomn:iomx));
            obs2(io) = sum( obs(1,iomn:iomx)...
                .* obs_wgt(iomn:iomx) ) ...
                / sum(obs_wgt(iomn:iomx));
        end
        DKL.(alg_name).obs{iNS} = obs2;
        DKL.(alg_name).obswt{iNS} = obs_wgt2;
        
        idx_below_05 = prob_bins < 0.5 & ~isnan(obs2);
        idx_above_05 = prob_bins >= 0.5 & ~isnan(obs2);
        obs_wgt2_temp = ones(size(obs_wgt2));
        DKL.(alg_name).PCAL2(iNS) = (sum(obs_wgt2_temp(idx_below_05).*(obs2(idx_below_05)-prob_bins(idx_below_05))) ...
                + sum(obs_wgt2_temp(idx_above_05).*(prob_bins(idx_above_05)-obs2(idx_above_05)))) ...
                / ( sum(obs_wgt2_temp(idx_below_05).*(0.5-prob_bins(idx_below_05)))...
                + sum(obs_wgt2_temp(idx_above_05).*(0.5-(1-prob_bins(idx_above_05)))));
        
        weight_tot = 0;
        N_times = 0;
        DKL.(alg_name).REL(iNS) = 0;
        for k1 = 1 : length(prob_bins)
            if ~isnan(obs2(k1))
                pk1 = prob_bins(k1);
                if 0.49 < pk1 && 0.51 > pk1
                    continue;
                end
                if obs2(k1)
                    N_times = N_times + 1;
                    DKL.(alg_name).REL(iNS) = ...
                        DKL.(alg_name).REL(iNS) ...
                        + obs_wgt2(k1)...
                        * KLD(pk1,obs2(k1));
                    weight_tot = weight_tot + obs_wgt2(k1);
                end
            end
        end
        DKL.(alg_name).REL(iNS) = DKL.(alg_name).REL(iNS) ...
            / weight_tot;    
            
    end
end

fig1 = figure(1);
metric_lists = {'correct','error','noclass','PCAL2'};
metric_axis_name_list = {'P_{Correct}','P_{Error}','P_{NoCall}','P_{Miscalibration}'};
metric_y_axis_list = {[0.5 1], [0 0.5], [0 0.5], [-1 1]};
metric_use_legend = [false true false false];
metric_use_exp = [false false false false];
metric_use_comp = [false false false false];
metric_y_tick_delta = [.1 .1 .1 .2];
alg_color = {[255 153 0]/255,'m',[0 204 255]/255,[0 128 0]/255};
alg_leg_name = {'DP ROC','HL NROC', 'BR', 'rBR'};
alg_2_compare = {'DP','HL_prob','BR','BRr'};
for im = 1 : length(metric_lists)
    subplot(2,2,im);
    hold off
    metric = metric_lists{im};
    metric_axis_name = metric_axis_name_list{im};
    metric_yaxis = metric_y_axis_list{im};
    plot([min(N_src_list) max(N_src_list)], [0 0],'k');
    hold on
    for ialg = 1 : length(alg_2_compare)
        alg_name = alg_2_compare{ialg};
        if metric_use_exp(im)
            h(ialg) = plot(N_src_list,exp(-DKL.(alg_name).(metric)),'color',alg_color{ialg},'linewidth',2);
        elseif metric_use_comp(im)
            h(ialg) = plot(N_src_list,1 - abs(DKL.(alg_name).(metric)),'color',alg_color{ialg},'linewidth',2);
        else
            h(ialg) = plot(N_src_list,DKL.(alg_name).(metric),'color',alg_color{ialg},'linewidth',2);
        end
    end
    grid on
    xlabel('Number of Sources');
    ylabel(metric_axis_name);
    if metric_use_legend(im)
        legend([h(:)], alg_leg_name,'location','best')
    end
    xticks([min(N_src_list):max(N_src_list)])
    yticks([min(metric_yaxis):metric_y_tick_delta(im):max(metric_yaxis)]);
    axis([min(N_src_list) max(N_src_list) min(metric_yaxis) max(metric_yaxis) ]);
end
saveas(fig1,'./plots/Pcorrect_Perror_Pnocall_Pcal.fig');
saveas(fig1,'./plots/Pcorrect_Perror_Pnocall_Pcal.png');

fig2 = figure(2);
alg_color = {[255 153 0]/255,'m',[0 204 255]/255,[0 128 0]/255};
alg_leg_name = {'DP ROC','HL NROC', 'BR', 'rBR'};
hold off
for ialg = 1 : length(alg_2_compare)
    alg_name = alg_2_compare{ialg};
    err_y = DKL.(alg_name).var_delta;
    h(ialg) = plot(N_src_list,DKL.(alg_name).avg_delta,'color',alg_color{ialg},'linewidth',2);
    hold on
    hr = plot([N_src_list N_src_list]',[-err_y err_y]'+[DKL.(alg_name).avg_delta DKL.(alg_name).avg_delta]','-r');
end
grid on
xlabel('Number of Sources');
ylabel({'max(Difference In Prediction)','Based On Fusion Order Changes'});
legend([h(:); hr(1)], [alg_leg_name '95% Conf'],'location','best')
ylim([0 1])
xticks([min(N_src_list):max(N_src_list)])
saveas(fig2,'./plots/Effects_of_fusion_order.fig');
saveas(fig2,'./plots/Effects_of_fusion_order.png');

fig3 = figure(3);
hold off
ialg_order = [4  3 1 2];
    
for ialg_temp = 1 : length(alg_2_compare)
    ialg = ialg_order(ialg_temp);
    alg_name = alg_2_compare{ialg};
    obs = DKL.(alg_name).obs{iNS};
    obs_wgt = DKL.(alg_name).obswt{iNS};
    if 1 == ialg_temp
        p = plot(prob_bins,prob_bins,'k','linewidth',2);
        hold on
    end
    h(ialg) = scatter(prob_bins(obs_wgt>0),obs(obs_wgt>0)...
        ,2000*obs_wgt(obs_wgt>0)/sum(obs_wgt(obs_wgt>0)),'filled'...
        ,'MarkerFaceColor',alg_color{ialg});
    
end
grid on
xlabel('Predicted Probability');
ylabel({'Observed Frequency'});
legend([p; h(:)], ['Perfectly Calibrated' alg_leg_name ],'location','best')
title(sprintf('%d Fused Sources',N_sources));
saveas(fig3,'./plots/obs_vs_pred.fig');
saveas(fig3,'./plots/obs_vs_pred.png');
end

function DKL = KLD(p,o)
DKL = o(0<o) .* log(o(0<o)./(p(0<o)*0.98 + .01));
end