classdef possibility_class < matlab.mixin.Copyable
    properties
        measure = [];  % possibility measure
    end
    methods
        function BBA = iCF(self)
            
            % get initial possiblity measure
            temp_measure = self.measure;
            
            % get number of classes
            N_classes = length(temp_measure);
            
            % sort possiblity measure
            [sort_measure,isort] = sort(temp_measure,'descend');
                        
            % initialize
            BBA(N_classes,1) = BBA_class;
            ib = 0;            
            init_focal = false(N_classes,1);
            
            % loop through number of classes
            for i = 1 : N_classes
                if i < length(temp_measure)
                    mass = sort_measure(i) - sort_measure(i+1);
                else
                    mass = sort_measure(i) + (1 - sort_measure(1));
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
        end
    end
end