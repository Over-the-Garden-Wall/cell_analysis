function strength = sa2strength(sa)
    strength = log(1+max(double(sa) - 50, 0));
end