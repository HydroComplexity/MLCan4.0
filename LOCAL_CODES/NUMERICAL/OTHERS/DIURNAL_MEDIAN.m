function [var_med, var_iqr] = DIURNAL_MEDIAN (hours, hours_all, var_in)

    for tt = 1:length(hours)

        inds = find(hours_all==hours(tt) & ~isnan(var_in));

        var_med(tt) = median(var_in(inds));
        var_iqr(tt) = iqr(var_in(inds));

    end