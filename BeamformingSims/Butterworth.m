
F_Sample = 16000;
F_Cutoff = 8000; 
order = 10;

Wn = F_Cutoff / (F_Sample / 2);

[b, a] = butter(order, Wn);

freqz(b, a, 1024, F_Sample);
