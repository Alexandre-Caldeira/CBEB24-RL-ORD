function state = discretize_val(val, min_val, max_val, num_states)
    
    state = nan(size(val));
    for idx = 1:numel(val)
        state(idx) = round(num_states * (val(idx) - min_val) / (max_val - min_val));
        if state(idx) >= num_states
            state(idx) = num_states;
        elseif state(idx) <= 0
            state(idx) = 1;
        end
    end
end
