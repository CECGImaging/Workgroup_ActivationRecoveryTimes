function [Q,T] = udq2quattrans(UDQ)
  % UDQ2QUATTRANS
  %
  % [Q,T] = udq2quattrans(UDQ)
  %
  % Convert a dual quaternion to rotation stored as a quaternion and a translation 
  %
  % Inputs:
  %  UDQ  list of rigid transformations for each rotation stored as dual
  %    quaternions
  %    2 by 4 #rotations
  % Output:
  %  Q  list of rotations stored as quaternions, one for each rotation
  %    #rotations by 4 (1,i,j,k) 
  %  T  list of translations stored as vectors, one for each rotation
  %    #rotations by 3
  %
  
  % Adapted from: http://isg.cs.tcd.ie/kavanl/dq/dqconv.c
  % regular quaternion (just copy the non-dual part):
  Q = UDQ(1,:,:)';
  % Translation vector
  T(:,1) = 2*( ...
    -UDQ(2,1,:).*UDQ(1,2,:) + ...
     UDQ(2,2,:).*UDQ(1,1,:) - ...
     UDQ(2,3,:).*UDQ(1,4,:) + ...
     UDQ(2,4,:).*UDQ(1,3,:));
  T(:,2) = 2.0*( ...
    -UDQ(2,1,:).*UDQ(1,3,:) + ...
     UDQ(2,2,:).*UDQ(1,4,:) + ...
     UDQ(2,3,:).*UDQ(1,1,:) - ...
     UDQ(2,4,:).*UDQ(1,2,:);
  T(:,3) = 2.0*( ...
    -UDQ(2,1,:).*UDQ(1,4,:) - ...
     UDQ(2,2,:).*UDQ(1,3,:) + ...
     UDQ(2,3,:).*UDQ(1,2,:) + ...
     UDQ(2,4,:).*UDQ(1,1,:);
end
