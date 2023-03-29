classdef object_class < matlab.mixin.Copyable
    properties
        mapped_info = info_class.empty;                 % mapped and calibrated information from source to object of interest
        calibrated_info = info_class.empty;              % any additional calibration to information as necessary
        confusion_matrix = confusion_matrix_class.empty; % class confusion matrix for additional calibration
    end
end