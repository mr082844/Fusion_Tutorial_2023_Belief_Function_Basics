classdef info_class < matlab.mixin.Copyable
    properties
        BBA = BBA_class.empty;                 % BBA information
        probability = probability_class.empty; % probability information
        possibility = possibility_class.empty; % possibility information
    end
end