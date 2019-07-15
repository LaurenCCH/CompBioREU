function[]=compare_plots_3D(data_nums, scale_small_probs, bw, p_k, p_c, c)

bw_matrix=zeros(length(bw),length(scale_small_probs),length(data_nums));
scale_small_matrix=zeros(length(bw),length(scale_small_probs),length(data_nums));
data_nums_matrix=zeros(length(bw),length(scale_small_probs),length(data_nums));



for i=1:length(bw)
    for j=1:length(scale_small_probs)
        for k=1:length(data_nums)
            bw_matrix(i,j,k)=bw(i);
            scale_small_matrix(i,j,k)=scale_small_probs(j);
            data_nums_matrix(i,j,k)=data_nums(k);
        end
    end
end

p_k_vector=reshape(p_k, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);
p_c_vector=reshape(p_c, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);
c_vector=reshape(c, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);

bw_vector=reshape(bw_matrix, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);
scale_small_vector=reshape(scale_small_matrix, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);
data_nums_vector=reshape(data_nums_matrix, [1,length(data_nums)*length(scale_small_probs)*length(bw)]);


p_k_scaled=(p_k_vector-(min(p_k_vector)))./((max(p_k_vector)-min(p_k_vector))*(length(p_k_vector)/(length(p_k_vector)-1)));
cm_pk=jet(length(p_k_vector));
color_index_pk=floor(p_k_scaled*length(cm_pk(:,1)))+1;
color_pk=cm_pk(color_index_pk,:);

p_c_scaled=(p_c_vector-(min(p_c_vector)))./((max(p_c_vector)-min(p_c_vector))*(length(p_c_vector)/(length(p_c_vector)-1)));
cm_pc=jet(length(p_c_vector));
color_index_pc=floor(p_c_scaled*length(cm_pc(:,1)))+1;
color_pc=cm_pc(color_index_pc,:);

c_scaled=(c_vector-(min(c_vector)))./((max(c_vector)-min(c_vector))*(length(c_vector)/(length(c_vector)-1)));
cm_c=jet(length(c_vector));
color_index_c=floor(c_scaled*length(cm_c(:,1)))+1;
color_c=cm_c(color_index_c,:);



scatter3(bw_vector,scale_small_vector,data_nums_vector, 20, color_pk);
colorbar
figure
scatter3(bw_vector, scale_small_vector,data_nums_vector, 20, color_pc);
colorbar
figure
scatter3(bw_vector, scale_small_vector,data_nums_vector, 20, color_c);
colorbar


end

