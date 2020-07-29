function y = convert(x,type)

if isa(x,type)
    y = x;
    return
end
if ~isa(x,'double')
    error('x is assumed to be a double')
end

switch type
    case 'mp'
        s = 10^(8-round(log10(ceil(max(x(:)))))+1);
        y = mp(round(s*x))/s;
    case 'sym'
        s = 10^(8-round(log10(ceil(max(x(:)))))+1);
        y = sym(round(s*x))/s;
    case 'single'
        y = single(x);
end