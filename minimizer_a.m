function out = minimizer_a(FD,alpha, T)
    out = FD .* (1-2*alpha*min(1,T./abs(FD)));
end

