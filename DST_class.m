classdef DST_class < handle
    properties
        source = source_class.empty;
        object_fusion_results = info_class.empty;
        FOD_align = FOD_mapping_class.empty;
        parameters = [];
    end
    methods
        function initialize(self,varargin)
            % read in variables
            optin = {[],[],[],[],[]};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [DST_param_file,FOD_param_file,src_names,N_src_tracks,N_obj_tracks] = optin{:};
            else
                return;
            end

            % build default signal struct
            if ~isfield(self.parameters,'file') ...
                    || ~strcmp(DST_param_file,self.parameters.file) ...
                    || ~isfield(self.parameters,'default_evidence')

                % get DST parameters
                DST_parameter_file_content  = read_parameters(DST_param_file);
                self.parameters = DST_parameter_file_content.parameters;
                self.parameters.file = DST_param_file;
                
                % get FOD parameters
                self.FOD_align = FOD_mapping_class;
                self.FOD_align.initialize(FOD_param_file);
                FOD_parameter_file_content  = read_parameters(FOD_param_file);
                FOD_set = unique(struct2cell(FOD_parameter_file_content.parameters.(self.parameters.integrator.FOD_set)));
                self.parameters.integrator.N_classes = length(FOD_set);
            end

            % build source data structure
            self.build_source_data_struct(src_names,N_src_tracks,N_obj_tracks)

            % build object_fusion_results data structure
            self.build_fusion_results_data_struct(N_obj_tracks)

            % initialize source structure
            self.source = self.builtinCopy(self.parameters.default_source);

            % initialize fusion results structure
            self.object_fusion_results = self.builtinCopy(self.parameters.default_fusion_results);
        end

        function process_source_data(self)
            self.condition_source_data();
            self.fuse_source_data();
        end

        function condition_source_data(self)

            % get number of sources
            N_sources = length(self.source);

            % remove prior
            for is = 1 : N_sources
                for it = 1 : length(self.source(is).track)
                    if isempty(self.source(is).track(it).given_info.probability.likelihood)
                        self.source(is).track(it).given_info.probability.remove_prior();
                    end
                end
            end

            % calibrate source data at their level (map to proper class set)
            for is = 1 : N_sources
                self.calibrate(self.source(is).track,is);
            end

            % map data to objects of interest
            self.map_2obj();

            % calibrate source data at object of interest level
            for is = 1 : N_sources
                self.calibrate(self.source(is).object,is);
            end
        end

        function calibrate(self,tracks_handle,is)
            % determine if source tracks and get confusion matrix
            switch class(tracks_handle)
                case 'track_class'
                    % construct confusion matrix
                    temp_CM_poss = [self.source(is).confusion_matrix.truthData.possibility];
                    CM_poss = [temp_CM_poss.measure];
                case 'object_class'
                    % do nothing at this point
                otherwise
                    error('MATLAB:input','class calibration not supported');
            end

            % get number of tracks to calibrate
            N_tracks = length(tracks_handle);

            % loop through tracks
            for it = 1 : N_tracks

                % get temp data structure to calibrate
                trk = tracks_handle(it);

                % determine class operation
                switch class(trk)
                    case 'track_class'
                        % initialize
                        starting_info = trk.given_info;
                        flg_FOD_align = true;

                    case 'object_class'
                        % initialize
                        starting_info = trk.mapped_info;
                        flg_FOD_align = false;

                        % construct confusion matrix
                        if ~isempty(trk.confusion_matrix)
                            temp_CM_poss = [trk.confusion_matrix.truthData.possibility];
                            CM_poss = [temp_CM_poss.measure];
                        else
                            CM_poss = [];
                        end
                    otherwise
                        error('MATLAB:input','class calibration not supported');
                end

                if isempty(CM_poss) %|| isequal(CM_poss,eye(size(CM_poss)))
                    % nothing needs to be done except maybe FOD mapping
                    if flg_FOD_align
                        starting_FOD = self.parameters.sources.(self.source(is).name).FOD_set;
                        ending_FOD = self.parameters.integrator.FOD_set;
                        trk.calibrated_info = self.FOD_align.mapFOD(starting_info,starting_FOD,ending_FOD);
                    else
                        trk.calibrated_info = starting_info;
                    end

                    % be sure to get these forms of information
                    if isempty(trk.calibrated_info.BBA)
                        trk.calibrated_info.BBA = trk.calibrated_info.probability.ipig(self.parameters.integrator.ipig_method);
                    end
                    if isempty(trk.calibrated_info.possibility.measure)
                        trk.calibrated_info.possibility = trk.calibrated_info.BBA.CF();
                    end
                else
                    % initialize
                    N_srcClasses = size(CM_poss,1);
                    calibrated_info = info_class;
                    calibrated_info.possibility = possibility_class;
                    calibrated_info.possibility.measure = ...
                        zeros(N_srcClasses,1);

                    % get BBA if doesn't already exist
                    if isempty(starting_info.BBA)
                        starting_info.BBA = starting_info.probability.ipig(self.parameters.integrator.ipig_method);
                    end

                    % calibrate
                    for ict = 1 : N_srcClasses
                        for ifcl = 1 : length(trk.given_info.BBA)
                            focal = starting_info.BBA(ifcl).focal;
                            possTgivC = 1 - prod(1 - CM_poss(ict,focal));
                            calibrated_info.possibility.measure(ict) = ...
                                calibrated_info.possibility.measure(ict) ...
                                + possTgivC * starting_info.BBA(ifcl).mass;
                        end
                    end

                    % get rest of forms of information
                    if isempty(calibrated_info.BBA)
                        calibrated_info.BBA = calibrated_info.possibility.iCF();
                    end
                    if isempty(calibrated_info.probability) ...
                            || isempty(calibrated_info.probability.posterior)
                        calibrated_info.probability = ...
                            calibrated_info.BBA.pig(...
                            starting_info.probability.prior);
                    end

                    % FOD align if necessary
                    if flg_FOD_align
                        starting_FOD = self.parameters.sources.(self.source(is).name).FOD_set;
                        ending_FOD = self.parameters.integrator.FOD_set;
                        trk.calibrated_info = self.FOD_align.mapFOD(calibrated_info,starting_FOD,ending_FOD);
                    else
                        trk.calibrated_info = calibrated_info;
                    end
                end
            end
        end

        function map_2obj(self)

            % get number of sources
            N_sources = length(self.source);

            % loop through number of sources
            for is = 1 : N_sources
                if isempty(self.source(is).object_association_mat) ...
                        || isequal(eye(size(self.source(is).object_association_mat))...
                        ,self.source(is).object_association_mat)
                    for it = 1 : length(self.source(is).track)
                        self.source(is).object(it).mapped_info = self.source(is).track(it).calibrated_info;
                    end
                else
                    [N_obj,N_trk] = size(self.source(is).object_association_mat);
                    if N_obj && N_trk
                        % get source track normalized likelihoods
                        sourceTracks_CALinfo = [self.source(is).track.calibrated_info];
                        sourceTracks_CALinfo_prob = [sourceTracks_CALinfo.probability];
                        sourceTracks_CALinfo_prob_like = [sourceTracks_CALinfo_prob.likelihood];

                        % create norm factor for AM
                        normFactor = sum(self.source(is).object_association_mat,2);
                        normFactorMat = repmat(normFactor,1,N_obj);

                        % avoid divide by zeros when normalizing
                        normFactorMat(normFactorMat==0) = 1;
                        AM_normalized = self.source(is).object_association_mat ...
                            ./ normFactorMat;

                        % only AM rows that sum to less than 1 get unknown applied
                        unknown_probability = ones(self.parameters.integrator.N_classes,1) ...
                            / self.parameters.integrator.N_classes;
                        unknownFactor = 1 - normFactor;
                        unknownFactor(unknownFactor < 0) = 0;

                        % map probabilities from source tracks to objects
                        object_CALinfo_prob_like = ([AM_normalized unknownFactor] ...
                            * [sourceTracks_CALinfo_prob_like unknown_probability]')';

                        % loop through objects
                        for io = 1 : N_obj
                            % store data
                            self.source(is).object(io).mapped_info.probability.initialize(object_CALinfo_prob_like(:,io));

                            % convert to other information domains
                            self.source(is).object(io).mapped_info.BBA = ...
                                self.source(is).object(io).mapped_info.probability.ipig('consonant');
                            self.source(is).object(io).mapped_info.possibility = ...
                                self.source(is).object(io).mapped_info.BBA.CF();
                        end
                    end
                end
            end
        end

        function fuse_source_data(self)

            % get number of sources
            N_sources = length(self.source);

            % initialize
            self.object_fusion_results = self.builtinCopy(self.parameters.default_fusion_results);
            init_info_2_fuse(N_sources) = info_class;

            % populate information to fuse
            % loop through objects of interest
            for io = 1 : length(self.object_fusion_results)
                % initialize
                info_2_fuse = init_info_2_fuse;

                for is = 1 : N_sources
                    info_2_fuse(is) = self.source(is).object(io).calibrated_info;
                end

                % fuse info
                self.object_fusion_results(io).BBA = self.NROC(info_2_fuse...
                    ,self.parameters.integrator.fusion_method);

                % convert to other information domains
                self.object_fusion_results(io).probability = ...
                    self.object_fusion_results(io).BBA.pig();
                self.object_fusion_results(io).probability = ...
                    self.object_fusion_results(io).probability.apply_prior();
                self.object_fusion_results(io).possibility = ...
                    self.object_fusion_results(io).BBA.CF();
            end
        end

        function fused_BBA = NROC(self,info_2_fuse,fusion_method)

            % get number of sources
            N_sources = length(info_2_fuse);

            % return if N_source < 2
            if 2 > N_sources
                fused_BBA = info_2_fuse.BBA;
                return;
            end

            % get dimention of fusion cube
            dim_list = NaN(1,N_sources);
            for is = 1 : N_sources
                dim_list(is) = length(info_2_fuse(is).BBA);
            end

            % initialize fused mass amd focals
            fused_BBA = BBA_class.empty;
            fused_mass = ones(dim_list);
            agg_source_focal = cell(dim_list);
            agg_source_mass = cell(dim_list);
            fused_root_conflict = cell(dim_list);
            N_el_dim = numel(fused_mass);

            % reshape to an 1xN array
            fused_mass = reshape(fused_mass,1,N_el_dim);
            agg_source_focal = reshape(agg_source_focal,1,N_el_dim);
            agg_source_mass = reshape(agg_source_mass,1,N_el_dim);
            fused_root_conflict = reshape(fused_root_conflict,1,N_el_dim);

            % fuse in sources
            for is = 1 : N_sources
                % initialize source mass and focals
                source_mass = ones(dim_list);
                source_focal = cell(dim_list);

                % create variable indexing (see: https://www.mathworks.com/matlabcentral/answers/362211-variable-indexing-for-n-dimension-data)
                var_idx = cell(size(dim_list));
                var_idx(:) = {':'};

                % populate source mass and focals
                for ib = 1 : dim_list(is)
                    var_idx{is} = ib;
                    source_mass(var_idx{:}) = info_2_fuse(is).BBA(ib).mass;
                    source_focal(var_idx{:}) = {{info_2_fuse(is).BBA(ib).focal}};
                end

                % reshape to an 1xN array
                source_mass = reshape(source_mass,1,N_el_dim);
                source_focal = reshape(source_focal,1,N_el_dim);

                % fuse mass
                fused_mass = fused_mass .* source_mass;
                if 1 == is
                    agg_source_focal = source_focal;
                    agg_source_mass = num2cell(source_mass);
                else
                    for ie = 1 : N_el_dim
                        agg_source_focal{ie} = [agg_source_focal{ie} source_focal{ie}];
                        agg_source_mass{ie}  = [agg_source_mass{ie} source_mass(ie)];
                    end
                end
            end

            % prune small mass
            if 1e4 < N_el_dim
                e_val = -8;
                flg_going_up = false;
                flg_going_down = false;
                for ieval = 1 : 6
                    thresh = eval(sprintf('1e%d',e_val));
                    if 0.99 < sum(fused_mass(thresh < fused_mass))
                        e_val = e_val + 1;
                        flg_going_up = true;
                        if flg_going_down
                            e_val = e_val - 1;
                            break;
                        end
                    else
                        e_val = e_val - 1;
                        flg_going_down = true;
                        if flg_going_up
                            break;
                        end
                    end
                end
                thresh = eval(sprintf('1e%d',e_val));
                idx_remove_small_mass = thresh > fused_mass;
                fused_mass(idx_remove_small_mass) = [];
                agg_source_focal(idx_remove_small_mass) = [];
                fused_root_conflict(idx_remove_small_mass) = [];
                N_el_dim = length(fused_mass);
            end

            % determine root conflict, support, and weight
            root_conflict_support = cell(1,N_el_dim);
            root_conflict_weight = cell(1,N_el_dim);

            % loop through fusion elements
            for ie = 1 : N_el_dim

                % determine root conflicts
                root_conflict = NaN(1,N_sources);
                for is = 1 : N_sources
                    root_conflict(is) = bi2de(reshape(agg_source_focal{ie}{is}...
                        ,1,numel(agg_source_focal{ie}{is})));
                end
                root_conflict = unique(root_conflict);
                N_roots = length(root_conflict);
                flg_is_root_conflict = false(1,N_roots);
                itter = 0;
                while ~all(flg_is_root_conflict) && 1 < N_roots && itter <= 100
                    itter = itter + 1;
                    roots_matrix = triu(bitand(root_conflict',root_conflict),1);
                    flg_is_root_conflict = ~any(roots_matrix,1) & ~any(roots_matrix,2)';
                    new_root_conflict = unique(roots_matrix(roots_matrix>0))';
                    root_conflict = [new_root_conflict root_conflict(flg_is_root_conflict)];
                    N_roots = length(root_conflict);
                end

                % determine root conflicts (slower but matches white paper)
                %                 root_conflict2 = NaN(1,N_sources);
                %                 for is = 1 : N_sources
                %                     root_conflict2(is) = bi2de(reshape(agg_source_focal{ie}{is}...
                %                         ,1,numel(agg_source_focal{ie}{is})));
                %                 end
                %                 root_conflict2 = unique(root_conflict2);
                %                 N_roots = length(root_conflict2);
                %                 i_list = 1 : N_roots;
                %                 New_root_conflict2 = NaN(1,factorial(N_roots));
                %                 inr2 = 0;
                %                 for k = i_list
                %                     i_comb = nchoosek(i_list,k);
                %                     for ic = 1 : size(i_comb,1)
                %                         other_i_list = i_list;
                %                         other_i_list(i_comb(ic,:)) = [];
                %                         this_i_comb = i_comb(ic,:);
                %                         int_this = self.parameters.integrator.N_classes^2-1;
                %                         for it = this_i_comb
                %                             int_this = bitand(int_this,root_conflict2(it));
                %                         end
                %                         if 0 < int_this
                %                             flg_root = true;
                %                             for io = other_i_list
                %                                 int_other = bitand(int_this,root_conflict2(io));
                %                                 if 0 < int_other
                %                                     flg_root = false;
                %                                     break;
                %                                 end
                %                             end
                %                         else
                %                             flg_root = false;
                %                         end
                %                         if flg_root
                %                             inr2 = inr2 + 1;
                %                             New_root_conflict2(inr2) = int_this;
                %                         end
                %                     end
                %                 end
                %                 root_conflict2 = sort(New_root_conflict2(~isnan(New_root_conflict2)));
                %                 N_roots = length(root_conflict2);

                bi_root_conflicts = de2bi([root_conflict 2^(self.parameters.integrator.N_classes-1)]);
                bi_root_conflicts(end,:) = [];
                if self.parameters.integrator.flg.DP
                    bi_root_conflicts = any(bi_root_conflicts,1);
                    N_roots = 1;
                end
                fused_root_conflict{ie} = bi_root_conflicts;
                N_root_conflicts = N_roots;
                switch fusion_method
                    case 'DROC'
                        root_conflict_weight{ie} = zeros(1,N_root_conflicts);
                    case 'HLNROC'
                        % determine support for root conflicts
                        root_conflict_support{ie} = zeros(1,N_root_conflicts);
                        if 1 < N_root_conflicts
                            for is = 1 : N_sources
                                src_support = zeros(1,N_root_conflicts);
                                for irc = 1 : N_root_conflicts
                                    if ~isempty(fused_root_conflict{ie}(irc,:))...
                                            && any(bitand(fused_root_conflict{ie}(irc,:)...
                                            ,reshape(agg_source_focal{ie}{is}...
                                            ,1,numel(agg_source_focal{ie}{is}))))
                                        src_support(irc) = 1;
                                    end
                                end
                                if ~all(src_support)
                                    root_conflict_support{ie} = root_conflict_support{ie} ...
                                        + src_support;
                                end
                            end
                        end
                        if self.parameters.integrator.probability_wgts ...
                                && ~all(0 == root_conflict_support{ie} - min(root_conflict_support{ie}))
                            support = root_conflict_support{ie};
                            weights = support.^support;
                            %                     [support_sorted,isort_supt] = sort(support,'descend');
                            %                     [~,i_unsort] = sort(isort_supt);
                            %                     support_primary_sorted = support_sorted / sum(support_sorted);
                            %                     weights = ones(1,N_root_conflicts); % initialize
                            %                     for ircp = 1 : N_root_conflicts-1
                            %                         support_primary_sorted_temp = support_primary_sorted;
                            %                         support_primary_sorted_temp(1:ircp) = geomean(support_primary_sorted_temp(1:ircp));
                            %                         power = support_sorted(ircp) - support_sorted(ircp+1);
                            %                         if power
                            %                             weights = weights.*support_primary_sorted_temp.^(power*2);
                            %                         end
                            %                     end
                            %                     root_conflict_support{ie} = weights(i_unsort) / sum(weights);
                            root_conflict_support{ie} = weights / sum(weights);
                        end

                        % determine weight for root conflicts
                        root_conflict_weight{ie} = zeros(1,N_root_conflicts);
                        for irc = 1 : N_root_conflicts
                            root_conflict_weight{ie}(irc) = ...
                                root_conflict_support{ie}(irc)...
                                * sum(fused_root_conflict{ie}(irc,:));
                        end
                    case 'HDNROC'
                        % determine support for root conflicts
                        root_conflict_support{ie} = zeros(1,N_root_conflicts);
                        if 1 < N_root_conflicts
                            N_soures_supporting = 0;
                            for is = 1 : N_sources
                                src_support = zeros(1,N_root_conflicts);
                                for irc = 1 : N_root_conflicts
                                    if ~isempty(fused_root_conflict{ie}(irc,:))...
                                            && any(bitand(fused_root_conflict{ie}(irc,:)...
                                            ,reshape(agg_source_focal{ie}{is}...
                                            ,1,numel(agg_source_focal{ie}{is}))))
                                        src_support(irc) = agg_source_mass{ie}(is);
                                    end
                                end
                                if ~all(src_support)
                                    root_conflict_support{ie} = root_conflict_support{ie} ...
                                        + src_support;
                                    N_soures_supporting = N_soures_supporting + 1;
                                end
                            end
                            if length(root_conflict_support{ie}) > 1
                                root_conflict_support{ie} = N_soures_supporting*root_conflict_support{ie} / sum(root_conflict_support{ie});
                            end
                        end

                        % determine weight for root conflicts
                        root_conflict_weight{ie} = zeros(1,N_root_conflicts);
                        for irc = 1 : N_root_conflicts
                            root_conflict_weight{ie}(irc) = ...
                                root_conflict_support{ie}(irc)...
                                * sum(fused_root_conflict{ie}(irc,:));
                        end
                    case 'HBNROC'
                        % determine support for root conflicts
                        root_conflict_support{ie} = zeros(1,N_root_conflicts);
                        if 1 < N_root_conflicts
                            for is = 1 : N_sources
                                src_support = zeros(1,N_root_conflicts);
                                for irc = 1 : N_root_conflicts
                                    if ~isempty(fused_root_conflict{ie}(irc,:))...
                                            && any(bitand(fused_root_conflict{ie}(irc,:)...
                                            ,reshape(agg_source_focal{ie}{is}...
                                            ,1,numel(agg_source_focal{ie}{is}))))
                                        src_support(irc) = 1;
                                    end
                                end
                                if ~all(src_support)
                                    root_conflict_support{ie} = root_conflict_support{ie} ...
                                        + src_support;
                                end
                            end
                            support = root_conflict_support{ie};
                            weights = binocdf(support,sum(support),1/length(support));
                            %                             weights = binocdf(support-1,sum(support)-1,1/length(support));
                            root_conflict_support{ie} = weights / sum(weights);
                        end

                        % determine weight for root conflicts
                        root_conflict_weight{ie} = zeros(1,N_root_conflicts);
                        for irc = 1 : N_root_conflicts
                            root_conflict_weight{ie}(irc) = ...
                                root_conflict_support{ie}(irc)...
                                * sum(fused_root_conflict{ie}(irc,:));
                        end
                end
                % sort focals and mass into fused results
                for irc = 1 : N_root_conflicts
                    flg_BBA_found = false;
                    N_bba = length(fused_BBA);
                    for ib = 1 : N_bba
                        if isequal(fused_BBA(ib).focal',fused_root_conflict{ie}(irc,:))
                            flg_BBA_found = true;
                            break;
                        end
                    end
                    if flg_BBA_found
                        if 1 < N_root_conflicts
                            if sum(root_conflict_weight{ie})
                                fused_BBA(ib).mass = fused_BBA(ib).mass ...
                                    + fused_mass(ie) * root_conflict_weight{ie}(irc) ...
                                    / sum(root_conflict_weight{ie});
                            end
                        else
                            fused_BBA(ib).mass = fused_BBA(ib).mass ...
                                + fused_mass(ie);
                        end
                    else
                        ib = N_bba + 1;
                        if 1 < N_root_conflicts
                            if sum(root_conflict_weight{ie})
                                fused_BBA(ib) = BBA_class;
                                fused_BBA(ib).mass = fused_mass(ie) ...
                                    * root_conflict_weight{ie}(irc) ...
                                    / sum(root_conflict_weight{ie});
                                fused_BBA(ib).focal = logical(fused_root_conflict{ie}(irc,:))';
                                fused_BBA(ib).cardonality = sum(fused_BBA(ib).focal);
                            end
                        else
                            fused_BBA(ib) = BBA_class;
                            fused_BBA(ib).mass = fused_mass(ie);
                            fused_BBA(ib).focal = logical(fused_root_conflict{ie}(irc,:))';
                            fused_BBA(ib).cardonality = sum(fused_BBA(ib).focal);
                        end
                    end
                end
            end

            % normalize BBA to unity
            sum_BBA = sum([fused_BBA.mass]);
            if 1 > sum_BBA
                N_bba = length(fused_BBA);
                for ib = 1 : N_bba
                    fused_BBA(ib).mass = fused_BBA(ib).mass / sum_BBA;
                end
            end

            % sort in ascending cardonality
            [~,isort] = sort([fused_BBA(:).cardonality]);
            fused_BBA = fused_BBA(isort);
        end

        function fused_BBA = HLNROC(self,info_2_fuse)

            % map legacy code inputs to new code inputs
            fused_BBA = NROC(self,info_2_fuse,'HLNROC');
        end

        function build_source_data_struct(self,varargin)
            % read in variables
            optin = {[],[],[]};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [src_names,N_src_tracks,N_obj_tracks] = optin{:};
            else
                return;
            end

            % initialize sources
            N_sources = length(src_names);
            default_source(N_sources) = source_class;
            temp_object(N_obj_tracks) = object_class;
            for is = 1 : N_sources
                default_source(is).name = src_names{is};
                temp_track(N_src_tracks(is)) = track_class;
                default_source(is).track = temp_track;
                default_source(is).confusion_matrix = confusion_matrix_class;
                default_source(is).confusion_matrix.initialize(...
                    self.parameters.sources.(default_source(is).name).confusion_matrix...
                    ,self.parameters.sources.(default_source(is).name).trust_level);
                for its = 1 : N_src_tracks(is)
                    default_source(is).track(its).given_info = info_class;
                    default_source(is).track(its).given_info.BBA = BBA_class.empty;
                    default_source(is).track(its).given_info.probability = probability_class;
                    default_source(is).track(its).given_info.possibility = possibility_class;
                    default_source(is).track(its).calibrated_info = info_class;
                    default_source(is).track(its).calibrated_info.BBA = BBA_class.empty;
                    default_source(is).track(its).calibrated_info.probability = probability_class;
                    default_source(is).track(its).calibrated_info.possibility = possibility_class;
                end
                default_source(is).object = temp_object;
                for its = 1 : N_obj_tracks
                    default_source(is).object(its).mapped_info = info_class;
                    default_source(is).object(its).mapped_info.BBA = BBA_class.empty;
                    default_source(is).object(its).mapped_info.probability = probability_class;
                    default_source(is).object(its).mapped_info.possibility = possibility_class;
                    default_source(is).object(its).calibrated_info = info_class;
                    default_source(is).object(its).calibrated_info.BBA = BBA_class.empty;
                    default_source(is).object(its).calibrated_info.probability = probability_class;
                    default_source(is).object(its).calibrated_info.possibility = possibility_class;
                end
                clear temp_track;
            end

            % store default source
            self.parameters.default_source = default_source;
        end

        function build_fusion_results_data_struct(self,varargin)
            % read in variables
            optin = {[]};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [N_obj_tracks] = optin{:};
            else
                return;
            end

            self.parameters.default_fusion_results(N_obj_tracks) = info_class;
        end

        function theCopy = builtinCopy(self,original)
            theCopy = copy(original);
            N_copys = length(theCopy);
            field_list = fields(theCopy);
            for ic = 1 : N_copys
                for ifld = 1 : length(field_list)
                    field_name = field_list{ifld};
                    class_type = class(theCopy(ic).(field_name));
                    if contains(class_type,'class')
                        theCopy(ic).(field_name) = self.builtinCopy(theCopy(ic).(field_name));
                    end
                end
            end
        end
    end
end