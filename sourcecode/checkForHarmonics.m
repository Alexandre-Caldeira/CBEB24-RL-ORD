function hasHarmonic = checkForHarmonics(a, b)
    hasHarmonic = false; % Initialize as false
    for i = 1:length(a)
        for j = 1:length(b)
            if mod(b(j), a(i)) == 0 % Check if b(j) is a multiple of a(i)
                hasHarmonic = [b(j), a(i)];
                return; % Exit the function as soon as a harmonic is found
            end
        end
    end
end
