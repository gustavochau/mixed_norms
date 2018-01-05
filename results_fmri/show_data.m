clc;
clear;
close all;

metodos{1} = 'sra';

metodos{2} = 'newton';
metodos{3} = 'fib';
metodos{4} = 'stef';

%  metodos{5} = 'fib_stef';
%  metodos{6} = 'fib_newt';
tau = '132.789';

figure
hold on
for ii=1:length(metodos)
    load(['mri_fromzero_' tau '_' char(metodos{ii}) '.mat'])
    plot(tiempos)
end
hold off
legend(metodos)
