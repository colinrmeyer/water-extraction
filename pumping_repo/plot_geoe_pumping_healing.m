% Plot effects of winter, pumping, and healing after extraction
close all
% load Models/Helheim_refined_interior_pump_heal.mat
load Models/Helheim_11points_pump_heal.mat
md1=md;

% % Single point
% xp=299807; yp=-2576740; % Helheim confluence *B*
% xp=274213; yp=-2563810; % Interior point for pumping test
% [a,pos] = min(sqrt((md.mesh.x-xp).^2+(md.mesh.y-yp).^2));
%
% for i=1:length(md1.results.TransientSolution)
%     N1(:,i)=md1.results.TransientSolution(i).EffectivePressure;
%     vel1(:,i)=md1.results.TransientSolution(i).Vel;
%     N1p(i)=md1.results.TransientSolution(i).EffectivePressure(pos);
%     vel1p(i)=md1.results.TransientSolution(i).Vel(pos);
% end

% % 11 points
xp=[299807 308516 293048 280612 294682 296884 290196 282422 286457 285247 302502];
yp=[-2576740 -2577550 -2564710 -2578780 -2566250 -2572380 -2578010 -2579160 -2558610 -2579600 -2577520];
for j=1:length(xp)
    [a,pos] = min(sqrt((md.mesh.x-xp(j)).^2+(md.mesh.y-yp(j)).^2));
    for i=1:length(md1.results.TransientSolution)
        N1(:,i)=md1.results.TransientSolution(i).EffectivePressure;
        vel1(:,i)=md1.results.TransientSolution(i).Vel;
        N1p(i)=md1.results.TransientSolution(i).EffectivePressure(pos);
        vel1p(i)=md1.results.TransientSolution(i).Vel(pos);
    end

    % t1=1/24:1/24:270;
    t1=0:1:270;

    % Plot 11 lines
    figure(1);
    subplot(2,1,1)
    set(gca,'FontSize',14);hold on
%     plot(t1,N1p./1e6,'linewidth',2);ylabel('N (MPa)')
plot(t1,(N1p-N1p(1))./1e6,'linewidth',2);ylabel('{\Delta}N (MPa)')
    axis([0 270 -1 40])
    subplot(2,1,2)
    set(gca,'FontSize',14);hold on
%     plot(t1,vel1p,'linewidth',2);ylabel('vel (m yr^{-1})')
plot(t1,vel1p-vel1p(1),'linewidth',2);ylabel('{\Delta}vel (m yr^{-1})')
    axis([0 270 -100 100])
    xlabel('Day')

end

% figure;
% subplot(2,1,1)
% set(gca,'FontSize',14);hold on
% plot(t1,mean(N1./1e6,1),'linewidth',2);ylabel('mean N (MPa)')
% axis([0 270 2.4 3.3])
% subplot(2,1,2)
% set(gca,'FontSize',14);hold on
% plot(t1,mean(vel1,1),'linewidth',2);ylabel('mean vel (m yr^{-1})')
% axis([0 270 2920 3000])
% xlabel('Day')

% figure;
% subplot(2,1,1)
% plot(t1,N1p,'linewidth',2);ylabel('N (Pa)')
% title('At pumping site')
% subplot(2,1,2)
% plot(t1,vel1p,'linewidth',2);ylabel('vel (m yr^{-1})')
% xlabel('Day')


