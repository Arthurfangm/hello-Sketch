%对所有数据集的图像进行特征提取，选择的特征是2angles_1ratio(两个角度特征一个线段%长度的
%组合特征)
run('C:\Users\方蒙\Documents\MATLAB\vlfeat-0.9.21\toolbox\vl_setup'); %启动vl-feat
%%
%分别为thetas特征矩阵、varphis特征矩阵和lratios特征矩阵规定统一的直方图化的边界
%19+19+11=49，对一条线段来说，有一个49*1的特征向量
thetas_edges = 10:20:350;
varphis_edges = 5:10:175;
l_ratios_edges = 0.05:0.1:0.95;

%%
%（边缘图像）数据集的地址
d = 'D:\sketch_based_image_retrieval\0605\query_demon';
cd(d);
images = dir(['*.' 'gif']);
number_of_images = size(images, 1); %number_of_images为文件夹中图片的数目

%%

processing_im = 0; %用于记录正在处理第几张图片
features = struct;
for im = images'
    imname = im.name;
    processing_im = processing_im+1;
    disp(im.name); %输出正在处理的图片的名字，可以查看程序运行的进度
    features(processing_im).name = imname;
    dd = [d, '\'];
    direction_of_im = [dd, imname];%获取正要处理的边缘图像的完整路径
    im_original = imread(direction_of_im);
    edgeim = edge(im_original,'canny', [0.1 0.2], 1);
    
    %%
    %将边缘图像线段化
    [edgelist, labeledgeim] = edgelink(edgeim, 10);
    tol = 2;
    seglist = lineseg(edgelist, tol);
    
    %%
    %获取线段的向量表示，中点、以及长度
    num_cells = size(seglist, 2);%记录有几个cell
    processing_cell = 1;%记录正在处理的第几个cell
    midpoint_vector = struct;%用于存放线段的向量表示、长度、中点坐标
    l_midpoint_vector = 0;%用于记录中点的个数
    
    while processing_cell <= num_cells
        precentcell = seglist{1, processing_cell};%获取当前要处理的cell
        num_endpoints = size(precentcell, 1);%表示当前的cell中有几个端点
        from_endpoint = 1;%用于记录处理到cell中的第几行
        while from_endpoint < num_endpoints
            to_endpoint = from_endpoint+1;
            l_midpoint_vector = l_midpoint_vector + 1;
            midpoint_vector(l_midpoint_vector).x = (precentcell(to_endpoint,2) - precentcell(from_endpoint, 2));
            midpoint_vector(l_midpoint_vector).y = (precentcell(to_endpoint,1) - precentcell(from_endpoint, 1));
            midpoint_vector(l_midpoint_vector).midx = (precentcell(from_endpoint,2) + precentcell(to_endpoint, 2))/2;
            midpoint_vector(l_midpoint_vector).midy = (precentcell(from_endpoint,1) + precentcell(to_endpoint, 1))/2;
            midpoint_vector(l_midpoint_vector).l = sqrt((precentcell(from_endpoint,1)-precentcell(to_endpoint,1))^2 + (precentcell(from_endpoint,2)-precentcell(to_endpoint,2))^2);
            from_endpoint = from_endpoint + 1;
        end
        processing_cell = processing_cell + 1;
    end
    
    a = zeros(2,1);
    b = zeros(2,1);
    thetas = zeros(l_midpoint_vector-1, l_midpoint_vector);%两条原线段的夹角0—360度
    varphis = zeros(l_midpoint_vector-1, l_midpoint_vector);%原线段与线段中点连线的夹角0-180度
    l_ratios = zeros(l_midpoint_vector-1, l_midpoint_vector);%线段长度的比值
    perim_segment_features = zeros(46, l_midpoint_vector);
    
    for i = 1 : l_midpoint_vector
        for j = 1 : l_midpoint_vector
            if(i==j) 
                continue;
            end
            a(1) = midpoint_vector(i).x;%a为基准线段的向量表示
            a(2) = midpoint_vector(i).y;
            
            %获取theta角
            b(1) = midpoint_vector(j).x;%b为比较线段的向量表示
            b(2) = midpoint_vector(j).y;
            cos_theta = dot(a,b)/midpoint_vector(i).l/midpoint_vector(j).l;
            theta = acosd(cos_theta);
            if (b(1)*a(2)-b(2)*a(1)<0)
                theta = 360 - theta;
            end
            
            %获取varphi角
            b(1) = midpoint_vector(j).midx - midpoint_vector(i).midx;
            b(2) = midpoint_vector(j).midy - midpoint_vector(i).midy;
            cos_varphi = dot(a,b)/midpoint_vector(i).l/sqrt(b(1)^2+b(2)^2);
            varphi = acosd(cos_varphi);
            
            %获取线段长度的比值
            l_ratio = midpoint_vector(i).l/midpoint_vector(j).l;
            if(l_ratio>1)
                l_ratio=1/l_ratio;
            end
            
            if(i<j)
                thetas(j-1,i) = theta;
                varphis(j-1,i) = varphi;
                l_ratios(j-1,i) = l_ratio;
            else
                thetas(j,i) = theta;
                varphis(j,i) = varphi;
                l_ratios(j,i) = l_ratio;
            end       
        end
        
    end
    
    h_thetas = hist(thetas, thetas_edges);
    h_varphis = hist(varphis, varphis_edges);
    h_l_ratios = hist(l_ratios, l_ratios_edges);
    
    for i = 1:l_midpoint_vector
        perim_segment_features(:,i) = [h_thetas(:,i); h_varphis(:,i); h_l_ratios(:, i)];
    end
end

%对草图的特征矩阵进行fisher编码，转化为向量
    %利用gmm_model文件夹中训练好的模型参数将每幅图像的表示从矩阵46*N维的矩阵转变为一维向量
    covariances = load('D:\sketch_based_image_retrieval\0601\gmm_model\covariances.mat');
    means = load('D:\sketch_based_image_retrieval\0601\gmm_model\means.mat');
    priors = load('D:\sketch_based_image_retrieval\0601\gmm_model\priors.mat');
    
    means = means.means;
    covariances = covariances.covariances;
    priors = priors.priors;
    
    %编码函数
    query_vector = vl_fisher(perim_segment_features, means, covariances, priors);
    
    distance_collection = struct;%用于存放自然图像的名称以及该自然图像与草图的距离
    
    %导入存有自然图像fisher编码的mat文件
    natural = load('D:\sketch_based_image_retrieval\0601\test_on_mepg7\fishered_mpeg7.mat');
    natural = natural.fishered_features;
    num_natural_im = size(natural, 2); %表示有几张真实图像
    for i = 1 : num_natural_im
        distance_collection(i).name = natural(i).name;
        distance_1 = query_vector - natural(i).fisher_vector;
        distance_2 = distance_1 .^ 2;
        distance_3 = sum(distance_2);
        distance = sqrt(distance_3);
        distance_collection(i).dist = distance;
    end
    
    %%
    %根据distance_collection,对距离进行排序，获取距离最小的前10幅图像
    top_ten = struct;
    for i = 1 : 10%遍历10次，获得10幅距离最近的图像
        top_ten(i).dist = Inf;%先将准备要填入的位置的距离设为无穷大
        for j = 1 : num_natural_im
            if distance_collection(j).dist < top_ten(i).dist
                no_of_min = j; %记录最小的距离在distance_collection里的序号
                top_ten(i).name = distance_collection(no_of_min).name;
                top_ten(i).dist = distance_collection(no_of_min).dist;
            end
        end
        %%将这一轮找到的最小值的位置用无穷大替换，避免后续又找到这里
        distance_collection(no_of_min).dist = Inf;
    end

%%
%可视化显示
%load('fishered_features.mat');
%addpath ETHZ_renamed
direction = 'D:\sketch_based_image_retrieval\dataset\original';
cd(direction);
figure(1);
fig=1;
for i=1:10
    picname=top_ten(i).name;
    subplot(2,5,i);
    imshow(picname);
end

        
        