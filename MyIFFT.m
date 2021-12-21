function Temporel = MyIFFT(Frequentiel)	

% Calcul de la transformée inverse
	Temporel = ifft( fftshift(Frequentiel) ) * sqrt(length(Frequentiel));


