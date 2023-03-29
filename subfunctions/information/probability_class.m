classdef probability_class < matlab.mixin.Copyable
    properties
        posterior = [];  % postierior probability
        likelihood = []; % what's left after removing prior from posterior (normalized to unity)
        prior = [];      % prior used in postier calculation
    end
    methods
        function initialize(self,varargin)
            % read in variables
            optin = {[],[]};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [self.posterior,self.prior] = optin{:};
            else
                return;
            end
            
            % assume uniform if not provided
            if isempty(self.prior)
                self.prior = ones(size(self.posterior)) / size(self.posterior,1);
            end
            
            % add missing probability evenly over classes
            self.posterior = self.posterior + repmat(...
                (1 - sum(self.posterior,1))./size(self.posterior,1),size(self.posterior,1),1);
            
            % remove prior into likelihood
            self.remove_prior();
        end
        
        function BBA = ipig(self,varargin)
            % read in variables
            optin = {'consonant'}; % option is consonant or bayesian
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [method] = optin{:};
            else
                BBA = BBA_class.empty;
                return;
            end

            % get initial possiblity measure
            temp_likelihood = self.likelihood;

            % get number of classes
            N_classes = length(temp_likelihood);

            % initialize
            BBA(N_classes,1) = BBA_class;
            ib = 0;
            init_focal = false(N_classes,1);

            % inverse pignistic transform method
            switch method
                case 'consonant'
                    % sort possiblity measure
                    [sort_likelihood,isort] = sort(temp_likelihood,'descend');

                    % loop through N_classes
                    for i = 1 : N_classes
                        if i < length(temp_likelihood)
                            mass = (sort_likelihood(i) - sort_likelihood(i+1))*i;
                        else
                            mass = sort_likelihood(i)*i;
                        end
                        if 1e-8 < mass
                            ib = ib + 1;
                            BBA(ib).mass = mass;
                            BBA(ib).focal = init_focal;
                            BBA(ib).focal(isort(1:i)) = true;
                            BBA(ib).cardonality = i;
                        end
                    end

                    % remove unused elements
                    if N_classes ~= ib
                        BBA(ib+1:N_classes) = [];

                        % renormalize if necessary
                        sum_mass = sum([BBA.mass]);
                        if 1 > sum_mass
                            for ib = 1 : length(BBA)
                                BBA(ib).mass = BBA(ib).mass / sum_mass;
                            end
                        end
                    end
                case 'bayesian'
                    for ib = 1 : N_classes
                        mass = temp_likelihood(ib);
                        BBA(ib).mass = mass;
                        BBA(ib).focal = init_focal;
                        BBA(ib).focal(ib) = true;
                        BBA(ib).cardonality = 1;
                    end
                    
                    % renormalize if necessary
                    sum_mass = sum([BBA.mass]);
                    if 1 > sum_mass
                        for ib = 1 : length(BBA)
                            BBA(ib).mass = BBA(ib).mass / sum_mass;
                        end
                    end
            end
        end
        
        function remove_prior(self)
            if isempty(self.prior) || all(all(1e-8 > abs(self.prior * size(self.prior,2) - 1)))
                self.likelihood = self.posterior;
            else
                % get number of classes
                N_classes = length(self.prior);
                
                % initialize
                M = ones(N_classes,N_classes);
                
                % populate M
                for r = 1 : N_classes - 1
                    for c = 1 : N_classes
                        if r == c
                            M(r,c) = (self.posterior(r) - 1) * self.prior(c);
                        else
                            M(r,c) = self.posterior(r) * self.prior(c);
                        end
                    end
                end
                M_inv = M^-1;
                self.likelihood = M_inv(:,end);
            end
        end
        
        function self = apply_prior(self)
            if isempty(self.prior) || all(1e-8 > abs(self.prior - self.prior(1)))
                self.posterior = self.likelihood;
            else
                % apply Bayes Rule
                self.posterior = self.likelihood .* self.prior / ...
                    sum(self.likelihood .* self.prior);
            end
        end
    end
end