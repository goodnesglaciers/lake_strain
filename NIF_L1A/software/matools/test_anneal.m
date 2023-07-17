%% to examine simulated annealing results
clear all

bounds = [-2 2; -2 2];
OPTIONS = [];

[mhat,F,model,energy,count]=anneal('rosenbrock',bounds, OPTIONS);

[nt, nr] = size(energy);
np = size(model,1)

E = reshape(energy, nt*nr, 1);
mods = reshape(model, np, nt*nr);

I = E < 2*F;
figure
plot(mods(1,:), mods(2,:), 'o')
hold on
plot(mods(1,I), mods(2,I), 'r+')

mhat
