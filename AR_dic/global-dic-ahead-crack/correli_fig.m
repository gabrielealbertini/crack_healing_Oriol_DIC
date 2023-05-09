function correli_fig(Mesh,Udic,Rdic,m,n) 
%
mesh_plot(Mesh,'Deformed',Udic,'Field',Rdic,'Wireframe',false,'EdgeColor','none')
axis image ij
% axis([0 n 0 m])
%colorbar
xlabel('$x$ (pixels)')
ylabel('$y$ (pixels)')


end

