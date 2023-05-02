function sa = rectangular_array(J, dx, dy )
%RECTANGULAR_ARRAY(J, dx, dy) Generate a 2D array with four elements
%   Generate a rectangular array with:
%     - sensor spacing in the x direction of dx
%     - sensor spacing in the y direction of dy
%     - center of the array in the origin

% Calculate the sensor positions in a matrix with all x-axis locations in
% the first column and all y-axis locations in the second column.
% width = sqrt(J);
% 
% x_cor = (-(width-1)/2 : (width-1)/2)';
% y_cor = flip(x_cor);
% for i = 1:width
%     for j = 1:widht
%         
%     end
% end
if (floor(sqrt(J))~=sqrt(J))
   warning('Input J must be a square number, otherwise the array sensor would not be exact number of assign')
end
length=floor(sqrt(J)); p = zeros(length^2,2);
% fill p with coordinates:
for i=1:length
    for j=1:length
    p((i-1)*length+(j),:)= [dx*((1-length)/2+(i-1)), dy*((1-length)/2+(j-1))];
    end
end
sa = array(p, 'Rectangular array');
end

