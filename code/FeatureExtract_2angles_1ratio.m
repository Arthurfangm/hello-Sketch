%���������ݼ���ͼ�����������ȡ��ѡ���������2angles_1ratio(�����Ƕ�����һ���߶�%���ȵ�
%�������)
run('C:\Users\����\Documents\MATLAB\vlfeat-0.9.21\toolbox\vl_setup'); %����vl-feat
%%
%�ֱ�Ϊthetas��������varphis���������lratios��������涨ͳһ��ֱ��ͼ���ı߽�
%19+19+11=49����һ���߶���˵����һ��49*1����������
thetas_edges = 10:20:350;
varphis_edges = 5:10:175;
l_ratios_edges = 0.05:0.1:0.95;

%%
%����Եͼ�����ݼ��ĵ�ַ
d = 'D:\sketch_based_image_retrieval\0605\query_demon';
cd(d);
images = dir(['*.' 'gif']);
number_of_images = size(images, 1); %number_of_imagesΪ�ļ�����ͼƬ����Ŀ

%%

processing_im = 0; %���ڼ�¼���ڴ���ڼ���ͼƬ
features = struct;
for im = images'
    imname = im.name;
    processing_im = processing_im+1;
    disp(im.name); %������ڴ����ͼƬ�����֣����Բ鿴�������еĽ���
    features(processing_im).name = imname;
    dd = [d, '\'];
    direction_of_im = [dd, imname];%��ȡ��Ҫ����ı�Եͼ�������·��
    im_original = imread(direction_of_im);
    edgeim = edge(im_original,'canny', [0.1 0.2], 1);
    
    %%
    %����Եͼ���߶λ�
    [edgelist, labeledgeim] = edgelink(edgeim, 10);
    tol = 2;
    seglist = lineseg(edgelist, tol);
    
    %%
    %��ȡ�߶ε�������ʾ���е㡢�Լ�����
    num_cells = size(seglist, 2);%��¼�м���cell
    processing_cell = 1;%��¼���ڴ���ĵڼ���cell
    midpoint_vector = struct;%���ڴ���߶ε�������ʾ�����ȡ��е�����
    l_midpoint_vector = 0;%���ڼ�¼�е�ĸ���
    
    while processing_cell <= num_cells
        precentcell = seglist{1, processing_cell};%��ȡ��ǰҪ�����cell
        num_endpoints = size(precentcell, 1);%��ʾ��ǰ��cell���м����˵�
        from_endpoint = 1;%���ڼ�¼����cell�еĵڼ���
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
    thetas = zeros(l_midpoint_vector-1, l_midpoint_vector);%����ԭ�߶εļн�0��360��
    varphis = zeros(l_midpoint_vector-1, l_midpoint_vector);%ԭ�߶����߶��е����ߵļн�0-180��
    l_ratios = zeros(l_midpoint_vector-1, l_midpoint_vector);%�߶γ��ȵı�ֵ
    perim_segment_features = zeros(46, l_midpoint_vector);
    
    for i = 1 : l_midpoint_vector
        for j = 1 : l_midpoint_vector
            if(i==j) 
                continue;
            end
            a(1) = midpoint_vector(i).x;%aΪ��׼�߶ε�������ʾ
            a(2) = midpoint_vector(i).y;
            
            %��ȡtheta��
            b(1) = midpoint_vector(j).x;%bΪ�Ƚ��߶ε�������ʾ
            b(2) = midpoint_vector(j).y;
            cos_theta = dot(a,b)/midpoint_vector(i).l/midpoint_vector(j).l;
            theta = acosd(cos_theta);
            if (b(1)*a(2)-b(2)*a(1)<0)
                theta = 360 - theta;
            end
            
            %��ȡvarphi��
            b(1) = midpoint_vector(j).midx - midpoint_vector(i).midx;
            b(2) = midpoint_vector(j).midy - midpoint_vector(i).midy;
            cos_varphi = dot(a,b)/midpoint_vector(i).l/sqrt(b(1)^2+b(2)^2);
            varphi = acosd(cos_varphi);
            
            %��ȡ�߶γ��ȵı�ֵ
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

%�Բ�ͼ�������������fisher���룬ת��Ϊ����
    %����gmm_model�ļ�����ѵ���õ�ģ�Ͳ�����ÿ��ͼ��ı�ʾ�Ӿ���46*Nά�ľ���ת��Ϊһά����
    covariances = load('D:\sketch_based_image_retrieval\0601\gmm_model\covariances.mat');
    means = load('D:\sketch_based_image_retrieval\0601\gmm_model\means.mat');
    priors = load('D:\sketch_based_image_retrieval\0601\gmm_model\priors.mat');
    
    means = means.means;
    covariances = covariances.covariances;
    priors = priors.priors;
    
    %���뺯��
    query_vector = vl_fisher(perim_segment_features, means, covariances, priors);
    
    distance_collection = struct;%���ڴ����Ȼͼ��������Լ�����Ȼͼ�����ͼ�ľ���
    
    %���������Ȼͼ��fisher�����mat�ļ�
    natural = load('D:\sketch_based_image_retrieval\0601\test_on_mepg7\fishered_mpeg7.mat');
    natural = natural.fishered_features;
    num_natural_im = size(natural, 2); %��ʾ�м�����ʵͼ��
    for i = 1 : num_natural_im
        distance_collection(i).name = natural(i).name;
        distance_1 = query_vector - natural(i).fisher_vector;
        distance_2 = distance_1 .^ 2;
        distance_3 = sum(distance_2);
        distance = sqrt(distance_3);
        distance_collection(i).dist = distance;
    end
    
    %%
    %����distance_collection,�Ծ���������򣬻�ȡ������С��ǰ10��ͼ��
    top_ten = struct;
    for i = 1 : 10%����10�Σ����10�����������ͼ��
        top_ten(i).dist = Inf;%�Ƚ�׼��Ҫ�����λ�õľ�����Ϊ�����
        for j = 1 : num_natural_im
            if distance_collection(j).dist < top_ten(i).dist
                no_of_min = j; %��¼��С�ľ�����distance_collection������
                top_ten(i).name = distance_collection(no_of_min).name;
                top_ten(i).dist = distance_collection(no_of_min).dist;
            end
        end
        %%����һ���ҵ�����Сֵ��λ����������滻������������ҵ�����
        distance_collection(no_of_min).dist = Inf;
    end

%%
%���ӻ���ʾ
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

        
        