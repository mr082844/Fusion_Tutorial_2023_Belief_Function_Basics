function analyze_results(varargin)
% read in variables
optin = {'./tests/results_same_pmiscal.mat'};
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
for ifx = 1 : length(field_names)
    field_name = field_names{ifx};
    DKL.(field_name).DS = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).REL = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).RES = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).UNC = ones(length(N_src_list),length(Pmc_list))*-log(0.5);
    DKL.(field_name).confMat = cell(length(N_src_list),length(Pmc_list));
    DKL.(field_name).error = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).noclass = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).correct = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).PCAL = zeros(length(N_src_list),length(Pmc_list));
    DKL.(field_name).avg_delta = zeros(length(N_src_list),length(Pmc_list));
    for iNS = 1:length(N_src_list)
        N_sources = N_src_list(iNS);
        for iPmc = 1 : length(Pmc_list)
            Pmc = Pmc_list(iPmc);
            idx_rows = results.(field_name).Pmiscal == Pmc ...
                & results.(field_name).N_src == N_sources ...
                & ~isnan(results.(field_name).fused_prob(:,1))...
                & ~isnan(results.(field_name).fused_prob(:,2));
            N_rows.(field_name) = sum(idx_rows);
            
            prob = results.(field_name).fused_prob(idx_rows,1);
            pred = zeros(size(prob));
            pred(prob>0.51) = 1;
            pred(prob<0.49) = 2;
            DKL.(field_name).confMat{iNS,iPmc} = ...
                confusionmat(results.(field_name).truth(idx_rows)...
                ,pred);
            if ~isempty(DKL.(field_name).confMat{iNS,iPmc})
                DKL.(field_name).error(iNS,iPmc) = (DKL.(field_name).confMat{iNS,iPmc}(end-1,end) ...
                    + DKL.(field_name).confMat{iNS,iPmc}(end,end-1)) ...
                    / sum(sum(DKL.(field_name).confMat{iNS,iPmc}));
                DKL.(field_name).correct(iNS,iPmc) = (DKL.(field_name).confMat{iNS,iPmc}(end-1,end-1) ...
                    + DKL.(field_name).confMat{iNS,iPmc}(end,end)) ...
                    / sum(sum(DKL.(field_name).confMat{iNS,iPmc}));
                DKL.(field_name).noclass(iNS,iPmc) = 1 ...
                    - DKL.(field_name).correct(iNS,iPmc) ...
                    - DKL.(field_name).error(iNS,iPmc); 
            end
            for ir = find(idx_rows')
                truth_temp = false(1,2);
                truth_temp(results.(field_name).truth(ir)) = true;
                DKL.(field_name).DS(iNS,iPmc) = DKL.(field_name).DS(iNS,iPmc) ...
                    + KLD(results.(field_name).fused_prob(ir,:)...
                    ,truth_temp);
            end
            DKL.(field_name).DS(iNS,iPmc) = DKL.(field_name).DS(iNS,iPmc) ...
                / N_rows.(field_name);
            
            MC_pc_list = unique(results.(field_name).MC_pc(idx_rows,1));
            MC_t_list = unique(results.(field_name).MC_t(idx_rows,1));
            delta = NaN(length(MC_pc_list),length(MC_t_list),2);
            for iMC_pc = 1:length(MC_pc_list)
                MC_pc = MC_pc_list(iMC_pc);
                idx_MC_pc = results.(field_name).MC_pc == MC_pc;
                for iMC_t = 1 : length(MC_t_list)
                    MC_t = MC_t_list(iMC_t);
                    idx_MC_t = results.(field_name).MC_t == MC_t;
                    for ic = 1 : 2
                        idx_c = results.(field_name).truth == ic;
                        idx_test = idx_MC_pc & idx_MC_t & idx_rows & idx_c;
                        if any(idx_test)
                            delta(iMC_pc,iMC_t) = max(results.(field_name).fused_prob(idx_test,1)) ...
                                - min(results.(field_name).fused_prob(idx_test,1));
                        end
                    end
                end
            end
            DKL.(field_name).avg_delta(iNS,iPmc) = mean(mean(mean(delta,'omitnan'),'omitnan'),'omitnan');
            
            DKL.(field_name).ok_counts = zeros(2,length(0 : delta_p : 1));
            
            idx_rows = results.(field_name).Pmiscal == Pmc ...
                & results.(field_name).N_src == N_sources;
            N_rows.(field_name) = sum(idx_rows);
            
            for ir = find(idx_rows')
                for k1 = 1 : length(prob_bins)
                    k2 = length(prob_bins) - k1 + 1;
                    pk1 = prob_bins(k1);
                    if results.(field_name).fused_prob(ir,1) ...
                            > pk1 - delta_p/2 ...
                            && results.(field_name).fused_prob(ir,1) ...
                            < pk1 + delta_p/2
                        DKL.(field_name).ok_counts(results.(field_name).truth(ir),k1) = ...
                            DKL.(field_name).ok_counts(results.(field_name).truth(ir),k1) + 1;
                        break;
                    end
                end
            end
            DKL.(field_name).ok = zeros(size(DKL.(field_name).ok_counts));
            DKL.(field_name).ok = DKL.(field_name).ok_counts ...
                ./ repmat(sum(DKL.(field_name).ok_counts,1),2,1);
            DKL.(field_name).ok_weight = sum(DKL.(field_name).ok_counts,1) ...
                / sum(sum(DKL.(field_name).ok_counts));
            N_times = 0;
            weight_tot = 0;
            weight_area_tot = 0;
            for k1 = 1 : length(prob_bins)
                if ~isnan(DKL.(field_name).ok(1,k1))
                    pk1 = prob_bins(k1);
                    if 0.49 < pk1 && 0.51 > pk1
                        continue;
                    end
                    k2 = length(prob_bins) - k1 + 1;
                    pk2 = prob_bins(k2);
                    if DKL.(field_name).ok(1,k1)
                        N_times = N_times + 1;
                        DKL.(field_name).REL(iNS,iPmc) = ...
                            DKL.(field_name).REL(iNS,iPmc) ...
                            + DKL.(field_name).ok_weight(k1)...
                            * KLD(pk1,DKL.(field_name).ok(1,k1));
                        weight_tot = weight_tot + DKL.(field_name).ok_weight(k1);
                        weight_area_tot = weight_area_tot ...
                            + pk1*DKL.(field_name).ok_weight(k2);
                        if pk1 > 0.5
                            DKL.(field_name).PCAL(iNS,iPmc) = ...
                                DKL.(field_name).PCAL(iNS,iPmc) ...
                                + DKL.(field_name).ok_weight(k1) ...
                                * (pk1 - DKL.(field_name).ok(1,k1));
                        else
                            DKL.(field_name).PCAL(iNS,iPmc) = ...
                                DKL.(field_name).PCAL(iNS,iPmc) ...
                                + DKL.(field_name).ok_weight(k1) ...
                                * (DKL.(field_name).ok(1,k1) - pk1);
                        end
                    end
                    if ~isnan(DKL.(field_name).ok(2,k2)) && DKL.(field_name).ok(2,k2)
                        N_times = N_times + 1;
                         DKL.(field_name).REL(iNS,iPmc) = ...
                             DKL.(field_name).REL(iNS,iPmc) ...
                            + DKL.(field_name).ok_weight(k2)...
                            * KLD(pk2,DKL.(field_name).ok(2,k2));
                        weight_tot = weight_tot + DKL.(field_name).ok_weight(k2);
                        weight_area_tot = weight_area_tot ...
                            + pk2*DKL.(field_name).ok_weight(k1);
                        if pk2 > 0.5
                            DKL.(field_name).PCAL(iNS,iPmc) = ...
                                DKL.(field_name).PCAL(iNS,iPmc) ...
                                + DKL.(field_name).ok_weight(k2) ...
                                * (pk2 - DKL.(field_name).ok(2,k2));
                        else
                            DKL.(field_name).PCAL(iNS,iPmc) = ...
                                DKL.(field_name).PCAL(iNS,iPmc) ...
                                + DKL.(field_name).ok_weight(k2) ...
                                * (DKL.(field_name).ok(1,k2) - pk2);
                        end
                    end
                end
            end
            DKL.(field_name).REL(iNS,iPmc) = DKL.(field_name).REL(iNS,iPmc) ...
                / weight_tot;
            DKL.(field_name).PCAL(iNS,iPmc) = DKL.(field_name).PCAL(iNS,iPmc) ...
                / weight_area_tot;
        end
    end
end

fig1 = figure(1);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,exp(-DKL.(field_name).REL(ihx,:)),'-','linewidth',2);
        hold on
    end
    grid minor
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('P_{REL}');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig1,'./plots/P_rel.fig');
saveas(fig1,'./plots/P_rel.png');

fig2 = figure(2);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,exp(-DKL.(field_name).DS(ihx,:)),'-','linewidth',2);
        hold on
    end
    grid minor
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('P_{KLD}');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig2,'./plots/P_KLD.fig');
saveas(fig2,'./plots/P_KLD.png');

fig3 = figure(3);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,DKL.(field_name).error(ihx,:),'-','linewidth',2);
        hold on
    end
    grid minor
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('P_{Error}');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig3,'./plots/P_error.fig');
saveas(fig3,'./plots/P_error.png');

fig4 = figure(4);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,DKL.(field_name).correct(ihx,:),'-','linewidth',2);
        hold on
    end
    grid minor
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('P_{Correct}');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig4,'./plots/P_error.fig');
saveas(fig4,'./plots/P_error.png');

fig5 = figure(5);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,DKL.(field_name).noclass(ihx,:),'-','linewidth',2);
        hold on
    end
    grid minor
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('P_{NoCall}');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig5,'./plots/P_error.fig');
saveas(fig5,'./plots/P_error.png');

fig6 = figure(6);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,DKL.(field_name).PCAL(ihx,:),'-','linewidth',2);
        hold on
    end
    grid minor
    plot([-1 1],[0 0],'k-');
    plot([-1 1],[0 0],'k-');
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('Fused P_{miscal}');
    legend(field_names,'location','best')
    axis([-1 1 -1 1]);
end
saveas(fig6,'./plots/P_miscal.fig');
saveas(fig6,'./plots/P_miscal.png');

fig7 = figure(7);
for ihx = 1 : length(N_src_list)
    subplot(ceil(sqrt(length(N_src_list))),ceil(sqrt(length(N_src_list))),ihx)
    hold off
    for ifx = 1 : length(field_names)
        field_name = field_names{ifx};
        plot(Pmc_list,DKL.(field_name).avg_delta(ihx,:),'-','linewidth',2);
        hold on
    end
    grid minor
    plot([-1 1],[0 0],'k-');
    plot([-1 1],[0 0],'k-');
    title(sprintf('# Sources = %d',N_src_list(ihx)));
    xlabel('Expert P_{miscal}');
    ylabel('Avg Delta Order Change');
    legend(field_names,'location','best')
    axis([-1 1 0 1]);
end
saveas(fig7,'./plots/delta_prob.fig');
saveas(fig7,'./plots/delta_prob.png');
end

function DKL = KLD(p,o)
DKL = o(0<o) .* log(o(0<o)./(p(0<o)*0.98 + .01));
end