function [ax,t,T]=nestedTiledLayout(dimChildren,dimParent,ChildOpts,ParentOpts)
% Creates a tiled layout of tiled layouts.
% The nested layouts each have dimension [dimChildren(1), dimChildren(2)],
% [nrows, ncols], while the parent layout has dimension [dimParents(1), dimParents(2)]
%
% ChildOpts and ParentOpts are optional structures of parameters for the 
% respective tiled layouts.
%
% Based on this answer at MATLAB Central by Matt J (2023-11):
% https://uk.mathworks.com/matlabcentral/answers/2053617-common-legend-for-each-row-of-tiled-layout#answer_1362422
%
% Example:
% % 5 rows, each with a common legend in the middle
% close all
% % Make Layout Skeleton
% dimParent = [5,1];
% dimChildren = [1,4];
% ChildOpts = struct('TileSpacing','Compact');
% ParentOpts = struct('TileSpacing','Compact');
% [ax,t]=nestedTiledLayout(dimChildren,dimParent,ChildOpts,ParentOpts);
% 
% % Populate Axes
% for k=1:numel(ax); plot(ax(k), rand(20,2)); end 
% % Populate Legends
% for i=1:numel(t)
%        Strings=compose("sample %d.%d",i,1:2);
%        lh(i) = legend(t(i).Children(1),Strings,'Orientation','horizontal');
%        lh(i).Layout.Tile='north';
% end

if nargin < 4
    ParentOpts = struct([]);
    if nargin < 3
        ChildOpts = struct([]);
    end
end

    T = tiledlayout(dimParent(1),dimParent(2));
    set(T,ParentOpts);
    k=0; 
    for i=1:prod(dimParent)
       t(i)=tiledlayout(T,dimChildren(1),dimChildren(2));
       t(i).Layout.Tile=i;   
       for j=1:prod(dimChildren)
           k=k+1;
           ax(k)=nexttile(t(i));       
       end
       set(t(i),ChildOpts)
    end 
end