function grid=uniqueGrid(values,bounds)
%UNIQUEGRID Merge roundoff duplicates and preserve exact requested endpoints.
% Mixing log grids can produce 800 and 799.9999999999998; a cell that thin
% makes polynomial derivatives meaningless even though unique() keeps it.
wasRow=isrow(values);
grid=sort(min(max(double(values(:)),bounds(1)),bounds(2)));
tolerance=64*eps(max(abs(grid)));
grid=grid([true;diff(grid)>tolerance]);
grid([1,end])=bounds;
if wasRow,grid=grid.';end
end
