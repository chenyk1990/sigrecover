function [map]=sr_seis()

map = [[0.5*ones(1,40),linspace(0.5,1,88),linspace(1,0,88),zeros(1,40)]',[0.25*ones(1,40),linspace(0.25,1,88),linspace(1,0,88),zeros(1,40)]',[zeros(1,40),linspace(0.,1,88),linspace(1,0,88),zeros(1,40)]'];

end