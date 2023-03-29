classdef confusion_matrix_class < matlab.mixin.Copyable
    properties
        truthData = struct('BBA',BBA_class.empty...   % BBA information
            ,'probability',probability_class.empty... % probability information
            ,'possibility',possibility_class.empty);  % possibility information
        trustFactor = [];% level of trust for source
    end
    methods
        function initialize(self,varargin)
            % read in variables
            optin = {[],0};
            if nargin <= length(optin)+1
                optin(1:nargin-1) = varargin(:);
                [posterior,self.trustFactor] = optin{:};
            else
                return;
            end
            
            % initialize N classes
            [N_truth_classes,N_call_classes] = size(posterior);
            
            % normalize posterior by calls
            sum_calls = sum(posterior,1);
            posterior(:,0<sum_calls) = posterior(:,0<sum_calls) ...
                ./ repmat(sum_calls(0<sum_calls),N_truth_classes,1);
            
            % initialize info classes
            self.truthData(1:N_call_classes) = struct('BBA',BBA_class.empty...   % BBA information
                ,'probability',probability_class.empty... % probability information
                ,'possibility',possibility_class.empty);  % possibility information
            
            % ideal CM
            ideal_CM = eye(N_call_classes);
            for icc = 1 : N_call_classes
                % initialize
                self.truthData(icc).BBA = BBA_class;
                self.truthData(icc).probability = probability_class;
                self.truthData(icc).possibility = possibility_class;
                
                % convert to BBA
                self.truthData(icc).probability.initialize(posterior(:,icc))
                self.truthData(icc).BBA = self.truthData(icc).probability.ipig();
                
                % convert to possibility
                temp_poss = self.truthData(icc).BBA.CF();
                
                % force diagonals to be one
                idx_diag = ideal_CM(:,icc) == 1;
                temp_poss.measure(idx_diag) = ideal_CM(idx_diag,icc);
                
                % linear combination of idea and actual based on trust
                self.truthData(icc).possibility.measure = self.trustFactor(icc) ...
                    * ideal_CM(:,icc) + (1-self.trustFactor(icc)) * temp_poss.measure;
            end
        end
    end
end