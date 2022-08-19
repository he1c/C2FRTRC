function x=noisemix(M,N,a,v1,v2,type)
%GAUSSMIX - Gaussian noise mixture generator
%   Usage: x = gaussmix(M,N,a,var1,var2);
%   M:  Number of samples per signal to be generated
%   N:  Number of signals to be generated
%   a:  Mixing factor per columns: 0<=a<1, fa(x)=(1-a)*f1(x)+a*f2(x) 
%   v1: Variance of the Gaussian nonimpulsive pdf f1
%   v2: Variance of the Gaussian impulsive pdf f2
%   x:  Output impulsive noise signal matrix or vector
%
%   EXAMPLE:
%   --------
%      v1=1;
%      v2=100;
%      a=0.05;
%      M=1e6;
%      N=2;
%      x=gaussmix(M,N,a,v1,v2);
%      var(x)
%      v1*(1-a)+v2*a
%      plot(x)
%
%  Copyright (C) 2004 Charalampos C. Tsimenidis
%  First created: Wed Aug 07 11:22:31 BST 2002
%  Last modified: Wed Jun  8 00:46:56 BST 2005
%  Revision: 0.8

%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


Ls=floor(a*M);
x=zeros(M,N);
for i=1:N
    [dummy, index]=sort(rand(M,1));
    if strcmp(type,'saltandpepper')
        x(index(1:Ls),i)=2*(binornd(1,0.5,Ls,1)-0.5);
    else
        x(index(1:Ls),i)=sqrt(v2)*randn(Ls,1);
    end
    if strcmp(type,'gaussian')
        x(index(Ls+1:end),i)=sqrt(v1)*randn(M-Ls,1);
    elseif  strcmp(type,'binary')
        x(index(Ls+1:end),i)=(sqrt(v1)*2*(binornd(1,0.5,M-Ls,1)-0.5));
    elseif  strcmp(type,'uniform')
        x(index(Ls+1:end),i)=unifrnd (-sqrt(v1),sqrt(v1),M-Ls,1);
    elseif  strcmp(type,'laplace')
        b=sqrt(v1)/sqrt(2);      %根据标准差求相应的b
        a=rand(M-Ls,1)-0.5;    %生成(-0.5,0.5)区间内均匀分布的随机数列 (一万个数的行向量);
        x(index(Ls+1:end),i)=0-b*sign(a).*log(1-2*abs(a));
    end
end    