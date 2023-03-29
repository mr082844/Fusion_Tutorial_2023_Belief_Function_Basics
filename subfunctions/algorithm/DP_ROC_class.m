classdef DP_ROC_class < handle
    properties
        source = info_class.empty;
        fused_results = info_class;
    end
    methods
        function fuse(self)
            
            % get number of sources
            N_sources = length(self.source);
            
            % initialize Fusion Results
            fused_BBA = BBA_class;
            
            % loop through sources
            for is = 1 : N_sources
                if 1 == is
                    fused_BBA = self.source(is).BBA;
                else
                    N_BBA_f = length(fused_BBA);
                    N_BBA_s = length(self.source(is).BBA);
                    
                    result_BBA = BBA_class.empty;
                    
                    for ibf = 1 : N_BBA_f
                        for ibs = 1 : N_BBA_s
                            
                            fused_mass = fused_BBA(ibf).mass * self.source(is).BBA(ibs).mass;
                            if 0 < fused_mass
                                fused_focal = bitand(fused_BBA(ibf).focal,self.source(is).BBA(ibs).focal);
                                if ~any(fused_focal)
                                    fused_focal = bitor(fused_BBA(ibf).focal,self.source(is).BBA(ibs).focal);
                                end
                                N_BBA_r = length(result_BBA);
                                flg_found_BBA_r = false;
                                for ibr = 1 : N_BBA_r
                                    if isequal(fused_focal,result_BBA(ibr).focal)
                                        flg_found_BBA_r = true;
                                        break;
                                    end
                                end
                                if flg_found_BBA_r
                                    result_BBA(ibr).mass = result_BBA(ibr).mass ...
                                        + fused_mass;
                                else
                                    ibr = N_BBA_r + 1;
                                    result_BBA(ibr) = BBA_class;
                                    result_BBA(ibr).mass = fused_mass;
                                    result_BBA(ibr).focal = fused_focal;
                                    result_BBA(ibr).cardonality = sum(result_BBA(ibr).focal);
                                end
                            end
                        end
                    end
                    self.fused_results.BBA = result_BBA;
                end
            end
            self.fused_results.probability = self.fused_results.BBA.pig();
        end
    end
end