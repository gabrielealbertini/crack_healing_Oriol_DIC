function correli_fig(Mesh,Udic,Rdic,m,n) 
%
mesh_plot(Mesh,'Deformed',Udic,'Field',Rdic)
axis image ij
% axis([0 n 0 m])
colorbar
xlabel('y, pixels')
ylabel('x, pixels')

