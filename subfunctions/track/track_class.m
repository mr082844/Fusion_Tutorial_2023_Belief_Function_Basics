classdef track_class < matlab.mixin.Copyable
    properties
        given_info = info_class.empty;                   % information as provided by source
        calibrated_info = info_class.empty;              % calibrated source information
    end
end