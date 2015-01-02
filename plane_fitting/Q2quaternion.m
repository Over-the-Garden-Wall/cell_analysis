function quat = Q2quaternion(Q)

% t = Qxx+Qyy+Qzz (trace of Q)
% r = sqrt(1+t)
% w = 0.5*r
% x = copysign(0.5*sqrt(1+Qxx-Qyy-Qzz), Qzy-Qyz)
% y = copysign(0.5*sqrt(1-Qxx+Qyy-Qzz), Qxz-Qzx)
% z = copysign(0.5*sqrt(1-Qxx-Qyy+Qzz), Qyx-Qxy)


t = Q(1,1)+Q(2,2)+Q(3,3);
r = sqrt(1+t);
w = r/2;
x = sqrt(1+Q(1,1)-Q(2,2)-Q(3,3))/2 * sign(Q(3,2)-Q(2,3));
y = sqrt(1-Q(1,1)+Q(2,2)-Q(3,3))/2 * sign(Q(1,3)-Q(3,1));
z = sqrt(1-Q(1,1)-Q(2,2)+Q(3,3))/2 * sign(Q(2,1)-Q(1,2));

quat = [w x y z];

    
            

end