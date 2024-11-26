function f=plot_match_iterate(cellPos)

% B = results{2,2}.fluorophorePosPx(:,1);
f=figure('Position', [100, 100, 600, 100])
hold on
for i=1:length(cellPos)
    A = cellPos{i};
    ax = plot(1:A(end)+100,i*ones(1,length(1:A(end)+100)),'red-','linewidth',1);
    ax.Color(4) = 0.2;
%     hold on
    plot(A,i*ones(1,length(A)),'black|','linewidth',1)
    % B = sort(Xn{end})';
%     ax2 = plot(1:B(end)+5,zeros(1,length(1:B(end)+5)),'red-','linewidth',4)
%     ax2.Color(4) = 0.2
%     hold on
%     plot(B,zeros(1,length(B)),'|','linewidth',1)
%     xlabel('Position (px)')
% legend({'Molecule','Optical map'},'Interpreter','latex','location','eastoutside')

end
% ylim([-1 length(cellPos)+1])
set(gca,'YTickLabel',[]);
end
