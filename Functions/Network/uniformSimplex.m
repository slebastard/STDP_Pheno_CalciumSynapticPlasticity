function samples = uniformSimplex(n_samples, dim)
    a = rand(dim-1,n_samples);
    a = sort(a,1);
    a = [zeros(1,n_samples); a; ones(1,n_samples)];
    samples = circshift(a, -1, 1) - a;
    samples = samples(1:end-1,:);
end

