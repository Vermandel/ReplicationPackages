function [rounded] = rounder(to_round,num_round)

rounded = ceil((to_round-0.4999999/(10^num_round))*(10^num_round))/(10^num_round);