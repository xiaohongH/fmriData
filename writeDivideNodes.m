function [ XindexTovnoList ] = writeDivideNodes( nodes1,nodes2,hemi,subName)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
addpath('/usr/local/freesurfer/matlab');
% 把2d平面的上的(x,y)坐标，转换成surface上的vertex;
% 找到patch_coor.x坐标对应的索引；
% 根据索引找到patch_coor.vno的值，即具体的vertex的值，

fileName = ['/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/',subName,'/surf/',hemi,'.full.flat.patch.3d']
patch_coor = read_patch(fileName);

%将vertex的number转换成成xy坐标
vertex_index = find(patch_coor.vno == nodes1);
x1 = patch_coor.x(:,vertex_index:vertex_index);
y1 = patch_coor.y(:,vertex_index:vertex_index);

vertex_index = find(patch_coor.vno == nodes2);
x2 = patch_coor.x(:,vertex_index:vertex_index);
y2 = patch_coor.y(:,vertex_index:vertex_index);

% totalDistance=sum(bsxfun(@minus,[x1,y1],[x2,y2]).^2,2); 

%求一元一次方程式
% 把线段分成divideNum--n等分，一共有divideNum+1个点的坐标；


totalDistance = sqrt(abs(x1-x2).^2 + abs(y1-y2).^2);
% % divideNum = fix(totalDistance/2);%根据距离划分多少个roi
divideNum = 20;%固定分为20个roi
ousideNum = 10;

partDistance = totalDistance/divideNum;
partX = (x1-x2)/divideNum;
partY = (y1-y2)/divideNum;

xlist = []
ylist = []



for i = -ousideNum:divideNum+ousideNum+1
    xlist = [xlist,x1-i*partX];
    ylist = [ylist,y1-i*partY];
end

% 求出所有的vertex
XindexTovnoList = [];
for i = 1:divideNum+2*ousideNum+1
    x = xlist(:,i:i);
    y = ylist(:,i:i);
    contrastList = abs(patch_coor.x -x);
    sortconList = sort(contrastList);
    kcontrastValue = sortconList(:,1:100);
    Xvalue = [];
    XToYvalue = [];
    Xindex = [];
    for j = kcontrastValue
        [row,Xcolumn]=find(contrastList==j);
        for k = length(Xcolumn)
            Xvalue = [Xvalue;patch_coor.x(:,Xcolumn(:,k:k):Xcolumn(:,k:k))];
            Xindex = [Xindex;Xcolumn(:,k:k)];
            XToYvalue =[XToYvalue;patch_coor.y(:,Xcolumn(:,k:k):Xcolumn(:,k:k))];
        end
    end
    %找到和y坐标最近的值
    contrastList = abs(XToYvalue - y);
    [YconListrow,Yconlistcolumn]=find(contrastList==min(min(contrastList)));
    coorColumn = Xindex(YconListrow:YconListrow,:); 
    % 二维坐标
    tkrx = patch_coor.x(:,coorColumn:coorColumn);
    tkry = patch_coor.y(:,coorColumn:coorColumn);

    % 二维坐标对应的vertex
    XindexTovnoList = [XindexTovnoList,patch_coor.vno(:,coorColumn:coorColumn)];
end

XindexTovnoList = XindexTovnoList'
combineAllNodes = []
%%
% surffn = '/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/surf/rh.pial.asc'
surffn = ['/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/',subName,'/surf/',hemi,'.pial.asc']

[v,f]=freesurfer_asc_load(surffn); % load surface 

nv=size(v,1); % number of vertices 
for idx=1:numel(XindexTovnoList)
    centerNode = XindexTovnoList(idx)
    aroundidxs=surfing_circleROI(v,f,centerNode,[2,5]); % indices of nodes in the disc 

    nodesFileName = ['newNodes_',hemi,'_SurfCoord.',string(centerNode),'.1D'];
    nodesFileName = join(nodesFileName,string(''));
    outNodesFilePath = ['/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/',char(subName),'/SUMA/',char(nodesFileName)];
%     outNodesFileName = join(outNodesFilePath,string(''));
    aroundidxs = aroundidxs'
%     save(outNodesFilePath,'aroundidxs','-ascii');
    combineAllNodes = [combineAllNodes;int64(aroundidxs)];
end
[combineAllNodes,is]=sort(combineAllNodes,'ascend');

%%
% combineAllNodes = combineAllNodes'
outFileName =['/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/',subName,'/SUMA/nodelocations_',hemi,'.1D']
combineAllNodesFileName = ['/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/',subName,'/SUMA/allnodelocations_',hemi,'.1D']

save(outFileName,'XindexTovnoList','-ascii');

% save(combineAllNodesFileName,'combineAllNodes','-ascii');
f = fopen(combineAllNodesFileName,'wt'); 
fprintf(f,'%d\n',combineAllNodes); 
fclose(f);
end

