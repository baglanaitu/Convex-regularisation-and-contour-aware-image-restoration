function Frequentiel = MyFFTRI(Temporel)	
	
% Calcul de la transformée
	Frequentiel = fftshift( fft(fftshift(Temporel)) );
