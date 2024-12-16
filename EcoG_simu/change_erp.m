function modified_data = change_erp(x)

% function to change ERP by randomly shifting time and amplitude.
latency_range = -30:30;
amplitude_range = 0.75:0.01:1.25;


data = permute(x, [3,1,2]);

[trials, channels, timepoints] = size(data);


modified_data = data;
for trial = 1:trials
    latency_shift = latency_range(randi(length(latency_range)));
    amplitude_scale = amplitude_range(randi(length(amplitude_range)));

    if latency_shift > 0
        modified_data(trial, :, :) = circshift(data(trial, :, :), latency_shift, 3);
        %modified_data(trial, :, 1:latency_shift) = [mean(modified_data(trial, :, 1:latency_shift),3) + randn(latency_shift,1)]';
    elseif latency_shift < 0
        shift_abs = abs(latency_shift);
        modified_data(trial, :, :) = circshift(data(trial, :, :), latency_shift, 3);
        %modified_data(trial, :, end-shift_abs+1:end) =  [mean(modified_data(trial, :, 1:latency_shift),3) + randn(shift_abs,1)]';
    end
    modified_data(trial, :, :) = modified_data(trial, :, :) * amplitude_scale;
end
modified_data = permute(modified_data, [2,3,1]);
