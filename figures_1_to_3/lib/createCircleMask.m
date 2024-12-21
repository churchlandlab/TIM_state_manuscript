function mask = createCircleMask(maskdims, center_x, center_y,radius)
[x, y] = meshgrid(1:maskdims(2), 1:maskdims(1));
mask = (x - center_x).^2 + (y - center_y).^2 <= radius^2;
end