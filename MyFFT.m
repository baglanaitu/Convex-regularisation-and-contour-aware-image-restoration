function Frequentiel = MyFFT(Temporel)	
	
% Calcul de la transformée
	Frequentiel = fftshift( fft(Temporel) ) / sqrt(length(Temporel)) ;

