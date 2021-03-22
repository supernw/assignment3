set(0, 'defaultFigureWindowStyle', 'docked')
clear
clf
close all

addpath('./code/');

Part1_Func(0.1);

%Part 2 and 3a
PART23_Func(0.8, 40e-9);

%Part 3b
wbs = linspace(35e-9, 45e-9, 5);
curr = zeros(1, length(wbs));
for i= 1:length(wbs)
    curr(i) = PART23_Func(0.8, wbs(i));
end

gaps = (100e-9/2-wbs)*2;

figure('Name', 'Gap')
plot(gaps, curr);
xlabel('Gap')
ylabel('Current')
title('Current vs. Varrying Gap')