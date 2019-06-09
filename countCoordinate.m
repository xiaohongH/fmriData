addpath('/usr/local/freesurfer/matlab');
% 把2d平面的上的(x,y)坐标，转换成surface上的vertex;
% 找到patch_coor.x坐标对应的索引；
% 根据索引找到patch_coor.vno的值，即具体的vertex的值，

patch_coor = read_patch(string('/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/surf/rh.full.flat.patch.3d'))
[sur_v,suf_f] = read_surf(string('/home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/surf/rh.white'))

%将vertex的number转换成成xy坐标
vertex_index = find(patch_coor.vno == 42980);
x1 = patch_coor.x(:,vertex_index:vertex_index);
y1 = patch_coor.y(:,vertex_index:vertex_index);

vertex_index = find(patch_coor.vno == 53144);
x2 = patch_coor.x(:,vertex_index:vertex_index);
y2 = patch_coor.y(:,vertex_index:vertex_index);

% totalDistance=sum(bsxfun(@minus,[x1,y1],[x2,y2]).^2,2); 

%求一元一次方程式
% 把线段分成divideNum--n等分，一共有divideNum+1个点的坐标；


totalDistance = sqrt(abs(x1-x2).^2 + abs(y1-y2).^2);
divideNum = fix(totalDistance/2);
ousideNum = 15;

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
    contrastList = abs(patch_coor.x - x);
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
save /home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/SUMA/nodelocations_rh.1D -ascii XindexTovnoList
%%
% 求出所有的freesurfer上的坐标以及所有的AFNI上的坐标

% fsRASToafniLIP = [];
% fsRASList = [];
% for m = 1:divideNum+2
%     XindexTovno = XindexTovnoList(:,m:m);
%     
%     vnoTotkRAS = sur_v(XindexTovno+1:XindexTovno+1,:);
% 
%     command = 'mri_info --vox2ras-tkr /home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/mri/orig.mgz';
%     [status,Torig]= system(command,'-echo');
%     command = 'mri_info --ras2vox-tkr /home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/mri/orig.mgz';
%     [status,ras2voxtkr]= system(command,'-echo');
% 
%     Torig = double(str2num(Torig))
%     ras2voxtkr = double(str2num(ras2voxtkr))
% 
%     tkRASTovolumeIndex = inv(Torig)*[vnoTotkRAS 1]';
%     volumeIndexTofsRAS = inv(ras2voxtkr)*tkRASTovolumeIndex;
%     fsRASList = [fsRASList;volumeIndexTofsRAS']
%     volumeIndexTofsRAS(1:2,:) = -volumeIndexTofsRAS(1:2,:);
%     
%     fsRASToafniLIP = [fsRASToafniLIP;volumeIndexTofsRAS'];
%     
% end

%%
% roinum = 1;
% for n = roinum:roinum
%     x = XindexTovnoList(1:1,roinum);
%     save /home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/SUMA/nodelocations.1D -ascii x
% end
%%
% for n = 6:6
%     x = fsRASList(n:n,1:3);
%     save /home/hh/study/python/code/ml/nipype_tutorial/bingfreesurfer/ab_data/tmp/edit.dat -ascii x
% end

% for n = 7:7
%     x = fsRASToafniLIP(n:n,1:3);
%     save /home/hh/study/graduate/fmridata/all_subj/surfaceProj/afni/ab/afniTofreesurferRoi/ROI_file.txt -ascii x
% end