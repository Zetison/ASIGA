function ss = mapleToMATLAB(s)

p = 0;
ss = [];
rootFuncDetected = false;
counter = 1;
n = numel(s);
i = 0;
while i < numel(s)
    i = i + 1;
    if i+6 <= n && strcmp(s(i:i+6), 'RootOf(')
        rootFuncDetected = true;
        ss = [ss, 'roots(['];
        rootFunc_p = p;
        p = p + 1;
        i = i + 6;
        counter = counter + 7;
    elseif s(i) == '('
        p = p + 1;
        ss = [ss, s(i)];
        counter = counter + 1;
    elseif s(i) == ')'
        p = p - 1;
        ss = [ss, s(i)];
        counter = counter + 1;
    elseif i+6 <= n && strcmp(s(i:i+6), '* cg ^ ')
        ss = [ss, ', '];
        i = i + 7;
        counter = counter + 2;
    elseif i+3 <= n && strcmp(s(i:i+3), '* cg')
        ss = [ss, ', '];
        i = i + 4;
        counter = counter + 2;
    else
        ss = [ss, s(i)];
        counter = counter + 1;
    end
    if rootFuncDetected && rootFunc_p == p
        ss = [ss, '])'];
        rootFuncDetected = false;
        counter = counter + 2;
    end
end
