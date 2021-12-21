clc, close all, clear all
load DataOne
%load DataTwo

% Parameters
mu     = 1;
N      = 300; % iterations
T      = 0.1;
alpha  = 0.1;
mu_bar = mu/alpha;

% Filters
f_h    = [0 0 0; 0 -1 1; 0 0 0];
f_v    = [0 0 0; 0 -1 0; 0 1 0];

F_H    = MyFFT2RI(f_h, 256); 
F_V    = MyFFT2RI(f_v, 256); 
norm_d = abs(F_H).^2 + abs(F_V).^2;

% Inputs
h0     = MyFFT2RI(IR, 256);
y0     = MyFFT2(Data); 
inv    = MyIFFT2(y0./h0);

y      = {};
y{1}   = y0;
et     = [];
j      = 0;

for k=1:N
    tic;    % regulizers
    FD_h     = MyIFFT2(F_H.*y{k});
    FD_v     = MyIFFT2(F_V.*y{k});
    
    % Updating of a
    a_h      = MyFFT2( minimizer_a(FD_h, T, alpha) );
    a_v      = MyFFT2( minimizer_a(FD_v, T, alpha) );
    
    % Updating of x 
    y{k+1}   = ( (abs(h0).^2 + mu_bar.*norm_d).^(-1) ) .* (h0'.*y0 + mu_bar*(conj(F_H).*a_h + conj(F_V).*a_v ));
    
    % Plot
    result   = MyIFFT2(y{k+1});
    imagesc(abs(result));
    colormap('gray');colorbar
    axis('square','off')
    title("Huber, mu = " + mu + ", T = " + T + ", K = " + k)
    drawnow

    %% Speed of convergence
    if rem(k, 25) == 0
        j = j+1;
        conv(j) = abs(sum(sum((y{k+1}-y{k}).^2))); % Euclidean (L2)
    end
end


% Analyzing in freq.domain
figure()
x_axis = linspace(-0.5,0.5, 256);
y_axis = linspace(-0.5,0.5, 256);
imagesc(x_axis, y_axis, log(abs(MyFFT2(result))));
colormap('gray'); colorbar
title("Huber, mu = "+mu + ", T = " + T)


conv'

