%get the average scross scan sessions
ave_data = zeros(20, 5, 50);
for x = 1:20
    for y = 1:5
        for z = 1:50
            ROI_data = roi_data(x, y, z, :);
            ave_data(x, y, z) = mean(ROI_data);
        end
    end
end

%pairwise - temporal isc
pairwise_temporal_ISC = zeros(20, 20, 5);
for z = 1:5
    for x = 1:20
        for y = 1:20
            data1 = reshape(ave_data(x, z, :), [1, 50]);
            data2 = reshape(ave_data(y, z, :), [1, 50]);
            pairwise_temporal_ISC(x, y, z) = corr(data1.', data2.');
        end
    end
end

%leave one out - temporal isc
loo_temporal_ISC = zeros(20, 5);
for x = 1:5
    new_loop = ave_data(:, x, :);
    for y = 1:20
        data = new_loop(1, 1, :);
        data1 = reshape(new_loop(1, 1, :), [1, 50]);
        new_loop = new_loop(2:20, 1, :);
        data2 = reshape(mean(new_loop), [1, 50]);
        loo_temporal_ISC(y, x) = corr(data1.', data2.');
        new_loop = vertcat(new_loop, data);
    end
end

%leave one out - dynamic isc
loo_dynamic_ISC = zeros(20, 5, 10);
for x = 1:5
    roi_loop = ave_data(:, x, :);
    i = 1;
    for y = 1:10
        sliced_data = roi_loop(:, 1, i:i + 4);
        i = i + 5;
        for z = 1:20
            data = sliced_data(1, 1, :);
            data1 = reshape(data, [1,5]);
            sliced_data = sliced_data(2:20, 1, :);
            data2 = reshape(mean(sliced_data), [1, 5]);
            loo_dynamic_ISC(z, x, y) = corr(data1.', data2.');
            sliced_data = vertcat(sliced_data, data);
        end
    end
end

%first subject dynamic isc
first_subject = zeros(10, 1);
i = 1;
for x = 1:10
    data = ave_data(1, 1, i:i + 4);
    first_subject(x, 1) = mean(data);
    i = i + 5;
end


%leave one out - spacial isc
loo_spatial_ISC = zeros(20, 1);
ave_roi = zeros(20, 5);
for y = 1:20
    for x = 1:5
        data = ave_data(y, x, :);
        ave_roi(y, x) = mean(data);
    end
end
loop_data = ave_roi;
for x = 1:20
    data1 = loop_data(1, :);
    loop_data = loop_data(2:20, :);
    data2 = mean(loop_data);
    loo_spatial_ISC(x, 1) = corr(data1.', data2.');
    loop_data = vertcat(loop_data, data1);
end

%intrasubject correlation
intrasubject_temporal_ISC = zeros(20, 5);
for x = 1:20
    for y = 1:5
        data1 = reshape(roi_data(x, y, :, 1), [1, 50]);
        data2 = reshape(roi_data(x, y, :, 2), [1, 50]);
        intrasubject_temporal_ISC(x, y) = corr(data1.', data2.');
    end
end

%inter-subject correlation
loo_ISFC = zeros(20, 5, 5);
for x = 1:5
    loop_data = ave_data;
    for y = 1:20
        data = loop_data(1, :, :);
        data1 = reshape(data(1, x, :), [1,50]);
        loop_data = loop_data(2:20, :, :);
        for z = 1:5 
            data2 = reshape(mean(loop_data(:, z, :)), [1,50]);
            loo_ISFC(y, x, z) = corr(data1.', data2.');
        loop_data = vertcat(loop_data, data);
        end
    end
end

%inter-subject representational similarity and behavior
%I calculate the similarity by first get the difference
%between two numbers. Then i normalize them into range 0 
%to 1 by using the formula (data - min(data)) / (max(data) - np.min(data))
%then I use 1 - the result of the above to get the similarity. 
behavior_similarity = zeros(20, 20);
for x = 1:20
    for y = 1:20
        behavior_similarity(x, y) = abs(behavior(1, x) - behavior(1, y));
    end
end
max_behav = max(behavior_similarity,[],'all');
for x = 1:20
    for y = 1:20
        behavior_similarity(x, y) = 1 - behavior_similarity(x, y)/ max_behav;
    end
end
%I hypothesize that the more similar between behaviors, 
%the stronger the intersubject correlation would be.
%So there is a correlation between behavior similarity
%and ISC. 
flatten_roi1 = untitled5(pairwise_temporal_ISC(:, :, 1));
flatten_roi2 = untitled5(pairwise_temporal_ISC(:, :, 2));
flatten_roi3 = untitled5(pairwise_temporal_ISC(:, :, 3));
flatten_roi4 = untitled5(pairwise_temporal_ISC(:, :, 4));
flatten_roi5 = untitled5(pairwise_temporal_ISC(:, :, 5));
flatten_behav = untitled5(behavior_similarity);
[rho1, pval1] = corr(flatten_roi1.', flatten_behav.','type','Spearman');
[rho2, pval2] = corr(flatten_roi2.', flatten_behav.','type','Spearman');
[rho3, pval3] = corr(flatten_roi3.', flatten_behav.','type','Spearman');
[rho4, pval4] = corr(flatten_roi4.', flatten_behav.','type','Spearman');
[rho5, pval5] = corr(flatten_roi5.', flatten_behav.','type','Spearman');
%based on the 5 p values, it seems that the ISC is not driven by the
%behavioral similarity as I have hypothesized. None of the 5 roi
%shows any statistically significance. 


    



            











