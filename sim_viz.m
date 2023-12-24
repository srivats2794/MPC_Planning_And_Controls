function sim_viz(pl,sim,viz,fbk,pl_rec) 
x_r_1 = [];
y_r_1 = [];

%% Reference Point and Pose
x1= sim.xf(1); y1= sim.xf(2); th1= sim.xf(3); 
R= [cos(th1), -sin(th1); sin(th1), cos(th1)]; 
yref_rect = [- viz.w/2,   viz.w/2, viz.w/2, - viz.w/2, - viz.w/2 ];
xref_rect = [- viz.l/2, - viz.l/2, viz.l/2,   viz.l/2, - viz.l/2];

T= R*[xref_rect;yref_rect];

x_f_rect= T(1,:)+x1; y_f_rect= T(2,:)+y1;

yref_tri = [0, -viz.marker_w/2, viz.marker_w/2, 0];
xref_tri = [viz.marker_h/2, - viz.marker_h/2, - viz.marker_h/2, viz.marker_h/2];
T= R*[xref_tri;yref_tri];

x_f_tri= T(1,:)+x1; y_f_tri= T(2,:)+y1;

f=figure(1);

set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);

for k = 1:length(fbk)
     
    subplot(1,2,1)
    plot(viz.map_x,viz.map_y,'-k','LineWidth',3);
    hold on
    fill(x_f_rect,y_f_rect,'k','HandleVisibility','off');
    hold on
    fill(x_f_tri,y_f_tri,'m','HandleVisibility','off');
    hold on

    
    for i=1:sim.obs_num
        fill(sim.obs_x(i)+viz.obs_circ_x, sim.obs_y(i)+viz.obs_circ_y,'r','HandleVisibility','off');
        hold on
    end
    
    x01 = fbk(1,k); y01 = fbk(2,k); th01 = fbk(3,k);

    if(norm(([x01;y01;th01]-sim.xf(1:3)),2) < 5e-3)
        break;
    end
    x_r_1 = [x_r_1 x01];
    y_r_1 = [y_r_1 y01];
    
    R= [cos(th01), -sin(th01); sin(th01), cos(th01)];
    T= R*[xref_rect;yref_rect];
    x0_rect= T(1,:)+x01; y0_rect= T(2,:)+y01;
    T= R*[xref_tri;yref_tri];
    x0_tri= T(1,:)+x01; y0_tri= T(2,:)+y01;
    
    r1=fill(x0_rect,y0_rect,'y','HandleVisibility','off');
    hold on
    r2=fill(x0_tri,y0_tri,'b','HandleVisibility','off');
    hold on
    plot(x01+viz.rob_circ_x,y01+viz.rob_circ_y,'--r','LineWidth',3);
    hold on
    r1.FaceColor='#7E2F8E';
    r1.EdgeColor='#7E2F8E';
    r2.FaceColor='#F56600';
    r2.EdgeColor='#F56600';
    
    plot(x_r_1,y_r_1,'Color','#7E2F8E','linewidth',3);hold on % plot exhibited trajectory
    hold on
    if k < length(fbk) % plot prediction
        plot( pl_rec(1:pl.N,1,k), pl_rec(1:pl.N,2,k),'--*','Color','#F56600')
    end  
    hold off
    %figure(500)
    title('Position Trace of Robot')
    ylabel('$Y$-position (m)','interpreter','latex','FontSize',20)
    xlabel('$X$-position (m)','interpreter','latex','FontSize',20)
    legend('Map Limits','Robot Obstacle Radius','Chosen Path','Predicted Path','FontSize',10, ...
        'Location','northwest')
    pause(0.1)
    
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end

close(gcf)

video = VideoWriter('exp.avi','Uncompressed AVI');

video = VideoWriter('exp.avi','Motion JPEG AVI');
video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
open(video)
writeVideo(video,F)
close (video)