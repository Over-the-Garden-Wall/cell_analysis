function Q = angles2Q(t, f, s)
    ct = cos(t); st = sin(t);
    cf = cos(f); sf = sin(f);
    cs = cos(s); ss = sin(s);

    Q = [ct*cs, sf*st*cs-cf*ss, sf*ss+cf*st*cs; ...
        ct*ss, cf*cs + sf*st*ss, cf*st*ss-sf*cs; ...
        -st, sf*ct, cf*ct];
end
    