close all
clear
clc

root_dir = fullfile('.','tests','results','calibration_study_con_50_amb_ipig_con_regression');
load(fullfile(root_dir,'results_test_IEEE_parametric_calibration_study_con_50_amb_ipig_con_regression.mat'));

flg_show_plots = false;
if flg_show_plots
    set(groot, 'DefaultFigureVisible', 'on');
else
    set(groot, 'DefaultFigureVisible', 'off');
end

%% set plots
flg_CM = 0;
flg_ROC = 0;
flg_pDist = 0;
flg_KLD = 1;
flg_wrapup = 1;
fusion_methods = {'HLnRoCprob','HLnRoCprob_DiscAsso','HLnRoCprob_noCal','HLnRoCprob_noCal_DiscAsso'...
    ,'HDnRoC','HDnRoC_DiscAsso','HDnRoC_noCal','HDnRoC_noCal_DiscAsso'...
    ,'HBnRoC','HBnRoC_DiscAsso','HBnRoC_noCal','HBnRoC_noCal_DiscAsso'...
    ,'DRoC','DRoC_DiscAsso','DRoC_noCal','DRoC_noCal_DiscAsso'};
nFusionMethods = length(fusion_methods);
source_names = {'Source1','Source2','Source3','Source4','Source5','Source6','Source7'};
N_sources = length(source_names);

%% create wrap up tables
if flg_wrapup
    tstart = tic();
    fprintf('creating wrapup tables...');

    % initialize wrap up table all data
    wrapup_filename = fullfile(root_dir,'wrap_up.xlsx');
    [N_tracks,N_classes,~,N_MC] = size(results.CPV_fuse);
    Alg_names = [fusion_methods source_names]';
    N_algs_w_sources = length(Alg_names);
    P_correct = NaN(N_algs_w_sources,1);
    P_error = NaN(N_algs_w_sources,1);
    P_NoCall = NaN(N_algs_w_sources,1);
    P_NoReport= NaN(N_algs_w_sources,1);
    P_CAL = NaN(N_algs_w_sources,N_classes);
    P_KLD = NaN(N_algs_w_sources,1);
    ROC_area = NaN(N_algs_w_sources,N_classes);
    wrap_up_table = table(Alg_names,P_correct,P_error,P_NoCall,P_NoReport,P_KLD,P_CAL,ROC_area);
    sheet_index = NaN(2,N_sources);
    sheet_index(1,:) = [2:2:N_sources*2+1];
    sheet_index(2,:) = [3:2:N_sources*2+1];
    
    % init for given/NOT sources
    P_correct = NaN(nFusionMethods+1,1);
    P_error = NaN(nFusionMethods+1,1);
    P_NoCall = NaN(nFusionMethods+1,1);
    P_NoReport= NaN(nFusionMethods+1,1);
    P_CAL = NaN(nFusionMethods+1,N_classes);
    P_KLD = NaN(nFusionMethods+1,1);
    ROC_area = NaN(nFusionMethods+1,N_classes);
    results_from = [fusion_methods,{'Source'}]';
    wrap_up_table_fuse_Src_int = table(results_from,P_correct,P_error,P_NoCall,P_NoReport,P_KLD,P_CAL,ROC_area);
    
    % populate using all data
    for ialg = 1 : nFusionMethods
        % Populate with algorithm information
        this_conf = results.conf_mat_fuse_agg(:,:,ialg);
        this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        wrap_up_table(ialg,2:end) = easyWrapUp(this_conf,this_CPV);
    end
    for is = 1 : N_sources
        % populate source in all data wrap up
        this_conf = results.conf_mat_src_agg(:,:,is);
        this_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        wrap_up_table(nFusionMethods+is,2:end) = easyWrapUp(this_conf,this_CPV);
    end
    sheet_name = sprintf('All_Data',is);
    writetable(wrap_up_table,wrapup_filename,'Sheet',sheet_name);
    
    % populate using source conditional data
    for is = 1 : N_sources
        % init fuse source related table
        wrap_up_table_given_src = wrap_up_table_fuse_Src_int;
        wrap_up_table_given_src(end,:) = wrap_up_table(nFusionMethods+is,:);
        wrap_up_table_given_NOT_src = wrap_up_table_fuse_Src_int;
        wrap_up_table_given_NOT_src(end,:) = wrap_up_table(nFusionMethods+is,:);
        for ialg = 1 : nFusionMethods
            % given source
            this_conf = results.conf_mat_fuse_src_agg(:,:,is,ialg);
            this_CPV = reshape(results.CPV_fuse_src(:,:,is,ialg,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
            wrap_up_table_given_src(ialg,2:end) = easyWrapUp(this_conf,this_CPV);
            
            % given not source
            this_conf = results.conf_mat_fuse_agg(:,:,ialg) - results.conf_mat_fuse_src_agg(:,:,is,ialg);
            this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
            source_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,~isnan(reshape(source_CPV(1,1,:),1,N_MC))) = [];
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,size(this_CPV,3)))) = [];
            wrap_up_table_given_NOT_src(ialg,2:end) = easyWrapUp(this_conf,this_CPV);
        end
        % save fused-source related tables
        sheet_name = sprintf('Given Source%d',is);
        writetable(wrap_up_table_given_src,wrapup_filename,'Sheet',sheet_name);
        sheet_name = sprintf('Given NOT Source%d',is);
        writetable(wrap_up_table_given_NOT_src,wrapup_filename,'Sheet',sheet_name);
    end
    tend = toc(tstart);
    fprintf(' complete (%0.2f sec)\n',tend);
end

%% Confusion matrix plots
if flg_CM
    tstart = tic();
    fprintf('Plotting Confusion Matrix...');
    for ialg = 1 : nFusionMethods
        % make save directory
        saveDir = fullfile(root_dir,'Confusion Matrix');
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% algorithm confusion matrix
        fig = figure('color','white');
        this_conf = results.conf_mat_fuse_agg(:,:,ialg);
        title_str = [fusion_methods{ialg},' Confusion Matrix'];
        easyPlotCM(this_conf,title_str);
        fig_name = fullfile(saveDir,[fusion_methods{ialg},' Confusion Matrix.png']);
        saveas(fig,fig_name);
        close(fig);
        
        for is = 1 : N_sources
            % make save directory
            saveDir = fullfile(root_dir,'Confusion Matrix',source_names{is});
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end

            %% algorithm confusion matrix given source reported
            fig = figure('color','white');
            this_conf = results.conf_mat_fuse_src_agg(:,:,is,ialg);
            title_str = {[fusion_methods{ialg},' Confusion Matrix'],['Given ', source_names{is},' Reported']};
            easyPlotCM(this_conf,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' Confusion Matrix Given ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
            %% algorithm confusion matrix given source didn't report
            fig = figure('color','white');
            this_conf = results.conf_mat_fuse_agg(:,:,ialg) - results.conf_mat_fuse_src_agg(:,:,is,ialg);
            title_str = {[fusion_methods{ialg},' Confusion Matrix'],['Given ', source_names{is},' Didn''t Report']};
            easyPlotCM(this_conf,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' Confusion Matrix Given Not ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
        end
    end
    
    for is = 1 : N_sources
        % make save directory
        saveDir = fullfile(root_dir,'Confusion Matrix',source_names{is});
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% Source confusion matrix
        fig = figure('color','white');
        this_conf = results.conf_mat_src_agg(:,:,is);
        title_str = [source_names{is},' Confusion Matrix'];
        easyPlotCM(this_conf,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Confusion Matrix'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source confusion matrix with discrete association
        fig = figure('color','white');
        this_conf = results.conf_mat_src_DA_agg(:,:,is);
        title_str = {[source_names{is},' Confusion Matrix'],'With Discrete Association'};
        easyPlotCM(this_conf,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Confusion Matrix DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source confusion matrix with no calibration
        fig = figure('color','white');
        this_conf = results.conf_mat_src_DA_agg(:,:,is);
        title_str = {[source_names{is},' Confusion Matrix'],'With No Calibration'};
        easyPlotCM(this_conf,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Confusion Matrix NoCal'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source confusion matrix with no calibration and discrete association
        fig = figure('color','white');
        this_conf = results.conf_mat_src_nCal_DA_agg(:,:,is);
        title_str = {[source_names{is},' Confusion Matrix'],'With No Calibration & Discrete Association'};
        easyPlotCM(this_conf,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Confusion Matrix NoCal and DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
    end
    close all
    tend = toc(tstart);
    fprintf(' complete (%0.2f sec)\n',tend);
end

%% ROC plots
if flg_ROC
    tstart = tic();
    fprintf('Plotting ROC...');

    [N_tracks,N_classes,~,N_MC] = size(results.CPV_fuse);
    for ialg = 1 : nFusionMethods
        % make save directory
        saveDir = fullfile(root_dir,'ROC Curve');
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% plot ROC curve for algorithm
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [fusion_methods{ialg},' Modified ROC Curve'];
        easyPlotMroc(this_CPV,title_str);
        fig_name = fullfile(saveDir,[fusion_methods{ialg},' Modified ROC Curve.png']);
        saveas(fig,fig_name);
        close(fig);
        
        for is = 1 : N_sources
            % make save directory
            saveDir = fullfile(root_dir,'ROC Curve',source_names{is});
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end

            %% algorithm ROC given source reported
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse_src(:,:,is,ialg,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
            title_str = {[fusion_methods{ialg},' Modified ROC Curve'],['Given ', source_names{is},' Reported']};
            easyPlotMroc(this_CPV,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' Modified ROC Curve Given ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
            %% algorithm ROC given source didn't report
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
            source_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,~isnan(reshape(source_CPV(1,1,:),1,N_MC))) = [];
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,size(this_CPV,3)))) = [];
            title_str = {[fusion_methods{ialg},' Modified ROC Curve'],['Given ', source_names{is},' Didn''t Report']};
            easyPlotMroc(this_CPV,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' Modified ROC Curve Given Not ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
        end
    end
    close all;
    
    for is = 1 : N_sources
        % make save directory
        saveDir = fullfile(root_dir,'ROC Curve',source_names{is});
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% Source ROC
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [source_names{is},' Modified ROC Curve'];
        easyPlotMroc(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Modified ROC Curve'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Modified ROC Curve'],'With Discrete Association'};
        easyPlotMroc(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Modified ROC Curve DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with no calibration
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Modified ROC Curve'],'With No Calibration'};
        easyPlotMroc(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Modified ROC Curve NoCal'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with no calibration and discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Modified ROC Curve'],'With No Calibration & Discrete Association'};
        easyPlotMroc(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' Modified ROC Curve NoCal and DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
    end
    tend = toc(tstart);
    fprintf(' complete (%0.2f sec)\n',tend);
end

%% probability distribution plots
if flg_pDist
    tstart = tic();
    fprintf('Plotting Probability Distribution...');

    [N_tracks,N_classes,~,N_MC] = size(results.CPV_fuse);
    for ialg = 1 : nFusionMethods
        % make save directory
        saveDir = fullfile(root_dir,'Pdistribution');
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% plot ROC curve for algorithm
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [fusion_methods{ialg},' Probability Distribution'];
        easyPlotpdist(this_CPV,title_str);
        fig_name = fullfile(saveDir,[fusion_methods{ialg},' pDist.png']);
        saveas(fig,fig_name);
        close(fig);
        
        for is = 1 : N_sources
            % make save directory
            saveDir = fullfile(root_dir,'Pdistribution',source_names{is});
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end

            %% algorithm ROC given source reported
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse_src(:,:,is,ialg,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
            title_str = {[fusion_methods{ialg},' Probability Distribution'],['Given ', source_names{is},' Reported']};
            easyPlotpdist(this_CPV,title_str);
            saveDir = fullfile(root_dir,'Pdistribution',source_names{is});
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' pDist Given ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
            %% algorithm ROC given source didn't report
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
            source_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,~isnan(reshape(source_CPV(1,1,:),1,N_MC))) = [];
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,size(this_CPV,3)))) = [];
            title_str = {[fusion_methods{ialg},' Probability Distribution'],['Given ', source_names{is},' Didn''t Report']};
            easyPlotpdist(this_CPV,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' pDist Given Not ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
        end
    end
    close all;
    
    for is = 1 : N_sources
        % make save directory
        saveDir = fullfile(root_dir,'Pdistribution',source_names{is});
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end
        
        %% Source ROC
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [source_names{is},' Probability Distribution'];
        easyPlotpdist(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' pDist'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Probability Distribution'],'With Discrete Association'};
        easyPlotpdist(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' pDist DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with no calibration
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Probability Distribution'],'With No Calibration'};
        easyPlotpdist(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' pDist NoCal'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source ROC with no calibration and discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Probability Distribution'],'With No Calibration & Discrete Association'};
        easyPlotpdist(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' pDist NoCal and DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
    end
    close all
    tend = toc(tstart);
    fprintf(' complete (%0.2f sec)\n',tend);
end

%% KLD plots
if flg_KLD
    tstart = tic();
    fprintf('Plotting KLD...');

    [N_tracks,N_classes,~,N_MC] = size(results.CPV_fuse);
    for ialg = 1 : nFusionMethods
        % make save directory
        saveDir = fullfile(root_dir,'KLD');
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% plot KLD curve for algorithm
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [fusion_methods{ialg},' Prediction vs Observation'];
        easyPlotKLD(this_CPV,title_str);
        fig_name = fullfile(saveDir,[fusion_methods{ialg},' KLD.png']);
        saveas(fig,fig_name);
        close(fig);
        
        for is = 1 : N_sources
            % make save directory
            saveDir = fullfile(root_dir,'KLD',source_names{is});
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end

            %% algorithm KLD given source reported
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse_src(:,:,is,ialg,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
            title_str = {[fusion_methods{ialg},' Prediction vs Observation'],['Given ', source_names{is},' Reported']};
            easyPlotKLD(this_CPV,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' KLD Given ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
            %% algorithm KLD given source didn't report
            fig = figure('color','white');
            this_CPV = reshape(results.CPV_fuse(:,:,ialg,:),N_tracks,N_classes,N_MC);
            source_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
            this_CPV(:,:,~isnan(reshape(source_CPV(1,1,:),1,N_MC))) = [];
            this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,size(this_CPV,3)))) = [];
            title_str = {[fusion_methods{ialg},' Prediction vs Observation'],['Given ', source_names{is},' Didn''t Report']};
            easyPlotKLD(this_CPV,title_str);
            fig_name = fullfile(saveDir,[fusion_methods{ialg},' KLD Given Not ',source_names{is},'.png']);
            saveas(fig,fig_name);
            close(fig);
            
        end
    end
    close all;
    
    for is = 1 : N_sources
        % make save directory
        saveDir = fullfile(root_dir,'KLD',source_names{is});
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end

        %% Source KLD
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src(:,:,is,:),N_tracks,N_classes,N_MC);
        this_CPV(:,:,isnan(reshape(this_CPV(1,1,:),1,N_MC))) = [];
        title_str = [source_names{is},' Prediction vs Observation'];
        easyPlotKLD(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' KLD'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source KLD with discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Prediction vs Observation'],'With Discrete Association'};
        easyPlotKLD(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' KLD DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source KLD with no calibration
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Prediction vs Observation'],'With No Calibration'};
        easyPlotKLD(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' KLD NoCal'],'.png']);
        saveas(fig,fig_name);
        close(fig);
        
        %% Source KLD with no calibration and discrete association
        fig = figure('color','white');
        this_CPV = reshape(results.CPV_src_nCal_DA(:,:,is,:),N_tracks,N_classes,N_MC);
        title_str = {[source_names{is},' Prediction vs Observation'],'With No Calibration & Discrete Association'};
        easyPlotKLD(this_CPV,title_str);
        fig_name = fullfile(saveDir,[[source_names{is},' KLD NoCal and DiscAsso'],'.png']);
        saveas(fig,fig_name);
        close(fig);
    end
    close all

    tend = toc(tstart);
    fprintf(' complete (%0.2f sec)\n',tend);
end
disp("done");

function wrap_up_table_entry = easyWrapUp(this_conf,this_CPV)
% initialize
[N_classes,~] = size(this_conf);
P_correct = NaN(1,1);
P_error = NaN(1,1);
P_NoCall = NaN(1,1);
P_NoReport= NaN(1,1);
P_CAL = NaN(1,N_classes);
P_KLD = NaN(1,1);
ROC_area = NaN(1,N_classes);
wrap_up_table_entry = table(P_correct,P_error,P_NoCall,P_NoReport,P_KLD,P_CAL,ROC_area);

% get basic information
N_correct = sum(diag(this_conf(1:N_classes,1:N_classes)));
N_total_calls = sum(sum(this_conf(1:N_classes,1:N_classes)));
N_NoCalls = sum(this_conf(:,N_classes+1));
N_NoReports = sum(this_conf(:,N_classes+2));
N_events = (N_total_calls+N_NoCalls+N_NoReports);

% Pcorrect
wrap_up_table_entry(1,1) = {N_correct / N_events};
% Perror
wrap_up_table_entry(1,2) = {(N_total_calls-N_correct) / N_events};
% P_NoCall
wrap_up_table_entry(1,3) = {N_NoCalls / N_events};
% P_NoReport
wrap_up_table_entry(1,4) = {N_NoReports / N_events};
% P_KLD
res = easyDetermineKLDmetrics(this_CPV);
wrap_up_table_entry(1,5) = {exp(-res.DS)};
% P_CAL
wrap_up_table_entry(1,6) = {exp(-max(0,res.REL)')};
% ROC area
[~,~,ROC_area] = easymROCcurve(this_CPV,1);
wrap_up_table_entry(1,7) = {ROC_area};
end

function easyPlotCM(this_confM,title_str)
[N_truth_classes,N_predicted_classes] = size(this_confM);
%     temp_conf = zeros(max(N_truth_classes,N_predicted_classes),max(N_truth_classes,N_predicted_classes));
%     temp_conf(1:N_truth_classes,1:N_predicted_classes) = this_conf;
%     this_conf = temp_conf;
total_entries = sum(sum(this_confM));
Y_truth = categorical(zeros(total_entries,1));
Y_pred = categorical(zeros(total_entries,1));
i = 0;
for itruth = 1 : N_truth_classes
    for ipred = 1 : N_predicted_classes
        i = i + 1;
        N_entries = this_confM(itruth,ipred);
        Y_truth(i:i+N_entries) = categorical(itruth-1);
        if ipred == N_truth_classes + 1
            Y_pred(i:i+N_entries) = categorical({'NoCall'});
        elseif ipred == N_truth_classes + 2
            Y_pred(i:i+N_entries) = categorical({'NoReport'});
        else
            Y_pred(i:i+N_entries) = categorical(ipred-1);
        end
        i = i + N_entries;
    end
end
plotconfusion(Y_truth,Y_pred);
title(title_str,'interp','none');
ylabel('Predicted Class','fontweight','bold','fontsize',12);
xlabel('Truth Class','fontweight','bold','fontsize',12);
end

function easyPlotMroc(this_CPV,title_str)
[~,N_classes,~] = size(this_CPV);
[fnr,fpr,~,~,~,truth,prediction] = easymROCcurve(this_CPV,0);
legend_str = [];
hold off;
plotroc(truth,prediction);
hold off;
colors = nncolor.ncolors(N_classes);
for ic = 1 : N_classes
    legend_str{end+1} = sprintf('Class %d',ic-1);
    plot(fnr{ic},fpr{ic},'LineWidth',2,'color',colors(ic,:))
    hold on;
end
line([1 0],[0 1],'LineWidth',2,'Color',[1 1 1]*0.8);
legend_str{end+1} = 'Random';

legend(legend_str,'location','best')
xlabel('False Negative Rate','fontweight','bold','fontsize',12);
ylabel('False Positive Rate','fontweight','bold','fontsize',12);
title(title_str);
grid on
end

function [fnr,fpr,area,tpr,thresholds,truth,prediction] = easymROCcurve(this_CPV,flg_area)
[N_tracks,N_classes,N_MC] = size(this_CPV);
prediction = zeros(N_classes,N_MC*N_tracks);
truth = zeros(N_classes,N_MC*N_tracks);
i = 0;
for it = 1 : N_tracks
    endi = i + N_MC;
    i = i + 1;
    track_CPV = reshape(this_CPV(it,:,:),N_classes,N_MC);
    iClass = mod(it-1,N_classes) + 1;
    truth(iClass,i:endi) = 1;
    prediction(:,i:endi) = track_CPV;
    i = endi;
end
[tpr,fpr,thresholds] = roc(truth,prediction);
fnr = cell(size(tpr));
for ic = 1 : N_classes
    fnr{ic} = 1-tpr{ic};
end
if ~isempty(flg_area) && flg_area
    area = zeros(1,N_classes);
    for ic = 1 : N_classes
        area(ic) = trapz(fpr{ic},tpr{ic});
    end
else
    area = NaN(1,N_classes);
end
end

function easyPlotpdist(this_CPV,title_str)
[N_tracks,N_classes,N_MCs] = size(this_CPV);
dim = ceil(sqrt(N_classes));
colors = nncolor.ncolors(N_classes);
prior_class_count = zeros(1,N_classes);
for it = 1 : N_tracks
    this_class = mod(it-1,N_classes)+1;
    prior_class_count(this_class) = prior_class_count(this_class)+1;
end
for ic = 1 : N_classes
    subplot(dim,dim,ic);
    hold off
    class_P = NaN(N_classes,N_MCs*max(prior_class_count));
    i = zeros(1,N_classes);
    for it = 1 : N_tracks
        this_class = mod(it-1,N_classes)+1;
        iend = i(this_class) + N_MCs;
        i(this_class) = i(this_class) + 1;
        class_P(this_class,i(this_class):iend) = this_CPV(it,ic,:);
        i(this_class) = iend;
    end
    legend_str = [];
    hold off
    for ic2 = 1 : N_classes
        [fx,x] = ecdf(class_P(ic2,:));
        plot([0; x; 1],[0; fx; 1],'color',colors(ic2,:),'linewidth',1);
        legend_str{end+1} = sprintf('Class %d', ic2 - 1);
        hold on
    end
    if ic == 4
        lh = legend(legend_str,'location','best','NumColumns',2);
        if isa(title_str,'cell')
            ntrows = length(title_str);
        else
            ntrows = 1;
        end
        if 1 == ntrows
            lh.Position(1) = 0.52 - lh.Position(3)/2;
            lh.Position(2) = 0.485 - lh.Position(4)/2;
        else
            lh.Position(1) = 0.52 - lh.Position(3)/2;
            lh.Position(2) = 0.458 - lh.Position(4)/2;
        end
    end
    xlabel({'Probability Given to',sprintf('Class %d',ic - 1)},'fontweight','bold','fontsize',10);
    ylabel('CDF','fontweight','bold','fontsize',12);
    yticks([0:.1:1])
    xticks([0:.1:1])
    xtickangle(90)
    sgtitle(title_str,'interp','none');
    grid on
end
end

function easyPlotKLD(this_CPV,title_str)
[~,N_classes,~] = size(this_CPV);
colors = nncolor.ncolors(N_classes);
res = easyDetermineKLDmetrics(this_CPV);
hold off
legend_str = [];
for ic = 1 : N_classes
    legend_str{end+1} = sprintf('Class %d (PCal = %0.2f)', ic-1,exp(-abs(res.REL(ic))));
    idx_not_zero = res.ok_weight(ic,:,ic)>0;
    scatter(res.prob_bins(idx_not_zero),res.ok(ic,idx_not_zero,ic)...
        ,max(5,(res.ok_weight(ic,idx_not_zero,ic)/sum(res.ok_weight(ic,idx_not_zero,ic))*125).^2),colors(ic,:),'filled'...
        ,'MarkerFaceAlpha',.6);
    hold on
end
plot([0 1],[0 1],'--k','linewidth',2);
legend_str{end+1} = 'Ideal';
axis([0 1 0 1])
xticks([0:.1:1]);
yticks([0:.1:1]);
legend(legend_str,'location','southoutside','NumColumns',3)
grid on
if isa(title_str,'cell')
    title_str = [title_str,{sprintf('Kullback–Leibler Divergence Probability is %0.2f',exp(-res.DS))}];
else
    title_str = [{title_str},{sprintf('Kullback–Leibler Divergence Probability is %0.2f',exp(-res.DS))}];
end
title(title_str)
xlabel('Predicted Probability','fontweight','bold','fontsize',12);
ylabel('Observed Frequency','fontweight','bold','fontsize',12);
end

function res = easyDetermineKLDmetrics(this_CPV)

res = [];
[N_tracks,N_classes,N_MC] = size(this_CPV);
prediction = zeros(N_classes,N_MC*N_tracks);
truth = false(N_classes,N_MC*N_tracks);

delta_p = 0.05;
prob_bins = 0 + delta_p/2 : delta_p : 1;
res.prob_bins = prob_bins;

i = 0;
for it = 1 : N_tracks
    endi = i + N_MC;
    i = i + 1;
    track_CPV = reshape(this_CPV(it,:,:),N_classes,N_MC);
    iClass = mod(it-1,N_classes) + 1;
    truth(iClass,i:endi) = true;
    prediction(:,i:endi) = track_CPV;
    i = endi;
end

% initialize
DS = 0;

[~,N_records] = size(truth);
for ir = 1 : N_records
    DS = DS + KLD(prediction(:,ir),truth(:,ir));
end
res.DS = DS / N_records;

ok_counts = zeros(N_classes,length(prob_bins),N_classes);
N_bins = length(prob_bins);
for ir = 1:N_records
    for ic = 1 : N_classes
        this_prob = prediction(ic,ir);
        k1 = this_prob >= prob_bins - delta_p/2 & this_prob < prob_bins + delta_p/2;
        ok_counts(truth(:,ir),k1,ic) = ok_counts(truth(:,ir),k1,ic) + 1;
    end
end
ok = ok_counts ./ repmat(sum(ok_counts,1),N_classes,1,1);
ok_weight = ok_counts ./ repmat(sum(sum(ok_counts,3),2),1,N_bins,N_classes);

res.ok = ok;
res.ok_weight = ok_weight;

REL = zeros(N_classes,1);

for ic = 1 : N_classes
    N_times = 0;
    weight_tot = 0;
    for k1 = 1 : length(prob_bins)
        if ~isnan(ok(1,k1,ic))
            pk1 = prob_bins(k1);
            if ok(ic,k1,ic)
                N_times = N_times + 1;
                REL(ic) = REL(ic) + sum(ok_weight(:,k1,ic)) * KLD(pk1,ok(ic,k1,ic));
                REL(ic) = REL(ic) + sum(ok_weight(:,k1,ic)) * KLD(1-pk1,1-ok(ic,k1,ic));
                weight_tot = weight_tot + sum(ok_weight(:,k1,ic));
            end
        end
    end
    REL(ic) = REL(ic) / weight_tot;
end
res.REL = REL;
end

function DKL = KLD(p,o)
if ~any(0<o)
    DKL = 0;
else
    DKL = o(0<o) .* log(o(0<o)./(p(0<o)*0.98 + 1/realmax));
end
end