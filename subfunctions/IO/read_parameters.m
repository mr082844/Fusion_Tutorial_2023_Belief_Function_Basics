function parameters = read_parameters(file_name)

% initialize
parameters = struct();

try
    fid = fopen(file_name);
    try
        i_line = 1;
        tline = fgetl(fid);
        nested_structs = cell(0);
        nested_structs_names = cell(0);
        layer_down = 0;
        while ischar(tline)
            try
                match = regexp(tline,'\s*(\S*)\s*(.*)%','tokens');
                if isempty(match)
                    match = regexp(tline,'\s*(\S*)\s*(.*)','tokens');
                end
                match = match{1};
                variable_name = match{1};
                variable_value = str2num(match{2});
                if ~isempty(match{2}) && isempty(variable_value)
                    str_list = strsplit(match{2},'''');
                    variable_value = str_list{2};
                end
                if isempty(variable_name)
                    % do nothing
                else
                    if isempty(variable_value) % structure or end
                        if strcmpi(variable_name,'end') && length(nested_structs) > 1
                            layer_down = layer_down - 1;
                            temp_struct = nested_structs{end-1};
                            temp_struct.(nested_structs_names{end}) = nested_structs{end};
                            nested_structs(end) = [];
                            nested_structs_names(end) = [];
                            nested_structs{end} = temp_struct;
                        elseif ~strcmpi(variable_name,'end')
                            layer_down = layer_down + 1;
                            nested_structs{end+1} = struct();
                            nested_structs_names{end+1} = variable_name;
                        else
                            layer_down = layer_down - 1;
                            parameters.(nested_structs_names{end}) = nested_structs{end};
                            nested_structs(end) = [];
                            nested_structs_names(end) = [];
                        end
                    else % variable
                        if ~isempty(nested_structs) && layer_down ~= 0
                            temp_struct = nested_structs{end};
                            temp_struct.(variable_name) = variable_value;
                            nested_structs{end} = temp_struct;
                        else
                            parameters.(variable_name) = variable_value;
                        end
                    end
                end
                i_line = i_line + 1;
                tline = fgetl(fid);
            catch
                warning('MATLAB:input','Error on line %d.',i_line);
            end
        end
        while layer_down > 0
            if strcmpi(variable_name,'end') && length(nested_structs) > 1
                layer_down = layer_down - 1;
                temp_struct = nested_structs{end-1};
                temp_struct.(nested_structs_names{end}) = nested_structs{end};
                nested_structs(end) = [];
                nested_structs_names(end) = [];
                nested_structs{end} = temp_struct;
            else
                layer_down = layer_down - 1;
                parameters.(nested_structs_names{end}) = nested_structs{end};
                nested_structs(end) = [];
                nested_structs_names(end) = [];
            end
        end
    catch 
        error('MATLAB:input','Could not fully read parameter file %s.',file_name);
    end
    fclose(fid);
catch
    error('MATLAB:input','Could not open parameter file %s.',file_name);
end
end