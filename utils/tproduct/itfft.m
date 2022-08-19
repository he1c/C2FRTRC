function [ t ] = itfft( t)
%ITFFT Summary of this function goes here
%   Detailed explanation goes here
% @param t: input/ouput data of tensor
% @date:19 Aug ,2018
s = size(t);
l=length(s);
for i = 3:l
    t=ifft(t,[],i);
end
end


