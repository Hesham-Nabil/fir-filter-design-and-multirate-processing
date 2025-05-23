% test_dft_comparison.m

N_values = [8, 16, 32, 64, 128, 256, 512, 1024];

for i = 1:length(N_values)
    N = N_values(i);
    fprintf('\n===== N = %d =====\n', N);
    compare_all(N);
end
