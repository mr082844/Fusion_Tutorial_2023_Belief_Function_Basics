function [x,y] = find_pCal_curve(pCal)
x = [0:0.01:1];
x(0.5 == x) = [];
if 0 == pCal
    y = x;
elseif 1 == abs(pCal)
    x = [0 x 1];
    x(x==0) = 1e-20;
    x(x==1) = 1 - 1e-20;
    y = x;
    y(1:length(y)/2) = 0.5-1e-9;
    y(length(y)/2+1:end) = 0.5+1e-9;
    y(1) = 1e-9;
    y(end) = 1-1e-9;
    if -1 == pCal
        x_temp = x;
        x = y;
        y = x_temp;
    end
else
    pCal_target = min(0.95,abs(pCal));
    persistent pCal_op_pt pCal_list;
    idx_op_pt = ismember(pCal_list,abs(pCal));
    if any(idx_op_pt)
        y = tanh(pCal_op_pt(idx_op_pt)*(x*2-1));
        y = y / max(y);
        y = (y + 1) / 2;
    else
        max_op = 1e5;
        min_op = 1e-20;
        op = mean([max_op,min_op]);
        pCal_list = [pCal_list abs(pCal)];
        y = tanh(op*(x*2-1));
        y = y / max(y);
        y = (y + 1) / 2;
%         plot(x,y);
        idx_below_05 = x < 0.5;
        idx_above_05 = x >=0.5;
        pCal_test = -(trapz(x(idx_below_05),y(idx_below_05)-x(idx_below_05)) ...
                + trapz(x(idx_above_05),x(idx_above_05)-y(idx_above_05))) ...
                / (trapz(x(idx_below_05),0.5 - x(idx_below_05))+trapz(x(idx_above_05),0.5-(1-x(idx_above_05))));
        if pCal_test < pCal_target
            min_op = op;
        else
            max_op = op;
        end
        op = mean([max_op,min_op]);
        while 1e-3 < abs(pCal_target - pCal_test)
            y = tanh(op*(x*2-1));
            y = y / max(y);
            y = (y + 1) / 2;
            %             plot(x,y);
            pCal_test = -(trapz(x(idx_below_05),y(idx_below_05)-x(idx_below_05)) ...
                + trapz(x(idx_above_05),x(idx_above_05)-y(idx_above_05))) ...
                / (trapz(x(idx_below_05),0.5 - x(idx_below_05))+trapz(x(idx_above_05),0.5-(1-x(idx_above_05))));
            if pCal_test < pCal_target
                min_op = op;
            else
                max_op = op;
            end
            op = mean([max_op,min_op]);
        end
        pCal_op_pt = [pCal_op_pt op];
    end
    if 0 < pCal
        temp_x = y;
        y = x;
        x = temp_x;
    end
end
end