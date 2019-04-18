function [p,q]=pade(f,xo,n,m,options,varargin)
% PADE computes the Padè approximant of the function F at the point XO.
%
% [P,Q]=PADE(F,XO,N,M) returns two polynomial forms, P of order N and Q of order M,
% representing the denominator and the numerator of the rational form 
% which approximates the function F around XO. F must accept X 
% and nth, F(X,nth). F return n-th derivative at X. 
% N and M must be positive integer values.
%
% [P,Q]=PADE(F,XO,N,M,OPTIONS) computed the Padè approximant with the 
% default optimization parameters replaced by values in the structure OPTIONS.
% Possible settings are: Vectorized{'off','on'}.
%
% [P,Q]=PADE(F,XO,N,M,OPTIONS,p1,p2,p3,...) computed the Padè approximant
% using p1, p2 and p3 as additional entries for F. 


% Reference paper: Baker, G. A. Jr. Essentials of Padé Approximants in 
% Theoretical Physics. New York: Academic Press, pp. 27-38, 1975.


% This routine has been programmed by Luigi Sanguigno, Phd.
% Affiliation: Italian Institute of Technology, Piazzale Tecchio
% 80, 80125, Naples, Italy. 
% First release 13-jun-2011. 

% EXAMPLE:
%    f=@(x,n) cat(2,log(1-x),-factorial(n(2:end)-1)./((1-x).^n(2:end)));
%    [p,q]=pade(f,0,6,6);
%    x=linspace(0,1,30);
%    plot(x,log(1-x),x,polyval(p,x)./polyval(q,x),'o');

fields={'Vectorized'};
values={{'on','off'}};

if (nargin<5)
    options=[];
end

% Assign default values to the fields not supplied by the user
for i=1:length(fields)
    if ~isfield(options,fields{i})
        options.(fields{i})=values{i}{1};
    end
end

% Calculate the n+m derivatives at xo
a=f(xo,0:(n+m)); 
a=(a./repmat(factorial(cast(0:(n+m),class(xo))),size(a,1),1)).';

% Calculation of Padè coefficients.
useHP = isa(xo,'mp');
if useHP
    d = mp.Digits;
else
    d = NaN;
end
pq = zeros(n+m+1,size(a,2),class(xo));
% for i = 1:size(a,2)
parfor i = 1:size(a,2)
    if useHP
        mp.Digits(d);
    end
    B = zeros(n+m+1,m,class(xo));
    for j = 1:m
        B(:,j) = [zeros(j,1); -a(1:n+m+1-j,i)];
    end
    pq(:,i) = cat(2,cat(1,speye(n+1),zeros(m,n+1)),B)\a(:,i);
%     pq(:,i) =cat(2,cat(1,speye(n+1),zeros(m,n+1)),...
%         spdiags(repmat(reshape(-cat(1,a(end:-1:1,i),0),1,[]),n+m+1,1),...
%         -(n+m+1):0,n+m+1,n))\a(:,i);
end

% Rewrite the output as evaluable polynomial forms.
p = zeros(n+1,size(a,2),class(xo));
q = zeros(m+1,size(a,2),class(xo));
% for i = 1:size(a,2)
parfor i = 1:size(a,2)
    if useHP
        mp.Digits(d);
    end
    p(:,i) = shiftpoly(pq(n+1:-1:1,i),xo); 
    q(:,i) = shiftpoly([pq(end:-1:n+2,i); ones(1,class(xo))],xo);
end


function ps = shiftpoly(p, xo)
% Displaces the origin of -xo

% Initialize values.
ps = zeros(size(p),class(xo)); 
q = ones(1,class(xo)); 
base = [1;-xo]; 
ps(end) = p(end);

% Substitute the base polynomial in the original polynomial form.
for n = 1:(length(p)-1)
    q = conv(q,base);
    ps(end-n:end) = ps(end-n:end)+p(end-n)*q;
end
    