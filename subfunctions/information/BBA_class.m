classdef BBA_class < matlab.mixin.Copyable
    properties
        mass = [];        % BBA mass
        focal = [];       % logical array, true for BBA focals (hyps with > 0 mass)
        cardonality = []; % number of elements in hyp
    end
    methods
        function probability_obj = pig(self,varargin)
            % read in variables
            optin = {[]};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [prior] = optin{:};
            else
                return;
            end
            
            % get number of tracks
            N_focals = length(self);
            
            % if no focals init to uniform probability
            if ~N_focals
                probability_obj = [];
            else
                % get number of classes                
                N_classes = length(self(1).focal);
                
                % initialize
                init_zeros = zeros(N_classes,1);
                probability_obj = probability_class;
                probability_obj.likelihood = init_zeros;
                
                % loop through focals
                for ib = 1 : N_focals
                    probability_obj.likelihood(self(ib).focal) = ...
                        probability_obj.likelihood(self(ib).focal) ...
                        + self(ib).mass / self(ib).cardonality;
                end

                % apply prior
                probability_obj.prior = prior;
                probability_obj.apply_prior();
            end            
        end
        
        function possibility_obj = CF(self)
            % get number of tracks
            N_focals = length(self);
            
            % if no focals init to uniform probability
            if ~N_focals
                possibility_obj = [];
            else
                % get number of classes                
                N_classes = length(self(1).focal);
                
                % initialize
                init_zeros = zeros(N_classes,1);
                possibility_obj = possibility_class;
                possibility_obj.measure = init_zeros;
                
                % loop through focals
                for ib = 1 : N_focals
                    possibility_obj.measure(self(ib).focal) = ...
                        possibility_obj.measure(self(ib).focal) ...
                        + self(ib).mass;
                end
            end
        end
    
        function self = consolidate(self)
            flg_consolidated = false(size(self));
            for ib = 2 : length(self)
                flg_BBA_found = false;
                for ibo = 1 : ib-1
                    if isequal(self(ibo).focal,self(ib).focal)
                        flg_BBA_found              = true;
                        break;
                    end
                end
                if flg_BBA_found
                    flg_consolidated(ib) = true;
                    self(ibo).mass       = self(ibo).mass + self(ib).mass;
                end
            end
            if any(flg_consolidated)
                self(flg_consolidated) = [];
            end
        end
    end
end