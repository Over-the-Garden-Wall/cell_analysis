function Q = quaternion2Q(quat)

w = quat(1);
x = quat(2);
y = quat(3);
z = quat(4);

Nq = w^2 + x^2 + y^2 + z^2;

if Nq > 0
    s = 2/Nq; 
else
    s = 0;
end

X = x*s; Y = y*s; Z = z*s;
wX = w*X; wY = w*Y; wZ = w*Z;
xX = x*X; xY = x*Y; xZ = x*Z;
yY = y*Y; yZ = y*Z; zZ = z*Z;

Q = ...
[ 1.0-(yY+zZ)       xY-wZ        xZ+wY  ; ...
      xY+wZ   1.0-(xX+zZ)       yZ-wX  ; ...
      xZ-wY        yZ+wX   1.0-(xX+yY) ];
end