% plot reachable set of y1-y2
fig = figure;
plot_2d_vertices 'build_48_y1_t.gen' 'b';
hold on; 
sn = @(x) 0.006;
fplot(sn,[0,20],'r');