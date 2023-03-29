classdef source_class < matlab.mixin.Copyable
    properties
        name = '';
        track                  = track_class.empty;            % tracked by source
        object                 = object_class.empty;           % object of interest (may not be source track)
        object_association_mat = [];                           % probability matrix relating source tracks to objects of interest
        confusion_matrix       = confusion_matrix_class.empty; % class confusion matrix for source
    end
end