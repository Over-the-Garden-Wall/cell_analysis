function Q = weighted_avg_Q(Q1,Q2,alpha)
    %average of 3d transforms Q1, Q2, with weight alpha 
    % Q = Q1^(alpha)*Q2^(1-alpha)
    
   %need to optimize, since the above computation will be slowish (I think)
   
   
  %find theta, phi, and psi for both Qs
  [theta1 phi1 psi1] = get_angles(Q1);
  [theta2 phi2 psi2] = get_angles(Q2);
  
  theta = theta1*alpha + theta2*(1-alpha);
  phi = phi1*alpha + phi2*(1-alpha);
  psi = psi1*alpha + psi2*(1-alpha);
  
  Q = angles2Q(theta, phi, psi);
  
end

    