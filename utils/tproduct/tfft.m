
function [t]= tfft(t)
%TFFT Summary of this function goes here
%   Detailed explanation goes here
% @param t: input/ouput data of tensor
% @date:19 Aug ,2018
% @author: haili
s = size(t);
l=length(s);
for i = 3:l
    t=fft(t,[],i);
end
end

