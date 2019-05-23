function luke_ProcessQPotLevelSets()
%save('GS_qpot_isosurfaces','fv1','fv2','fv3','fv4','fv5','fv6','fv7','fv8','fv9','fv10','fv11','isoval');
data = load('GS_qpot_isosurfaces.mat');
n0 = 500; % the desired number of faces
n1 = 2000;
figure;
hold on; grid;
view(3);
col = [winter(4);spring(4);[1 0 0];[0 1 0];[0 0 1]];

for k = 1 : 10

    %load up original faces	
	my_field =strcat('fv',num2str(k));
	s = data.(my_field);
	
	v = s.vertices;
    nv = size(v,1); 
    f = s.faces;
    nf = size(f,1); % the number of faces
    
    %smaller level sets need less points,
    %n0 < n1
    if ismember(k,[1,2,5,6])
        N = n0;
    else
        N = n1;
    end
        
    rfac = N/nf;
    [f0,v0] = reducepatch(f,v,rfac);
    fprintf('k = %d, nf = %d, nv = %d, nf0 = %d, nv0 = %d\n',k,nf,nv,size(f0,1),size(v0,1));
    p = patch('Vertices',v0,'Faces',f0,'Facecolor',col(k,:),'Edgecolor',col(k,:));
    alpha(0.3);
    
   

	%save reduced patches to output data
	my_field = strcat('f0_',num2str(k));
	reduced_patches.(my_field) = f0;
	my_field = strcat('v0_',num2str(k));
	reduced_patches.(my_field) = v0;

end

save('reduced_QPotLevelSets.mat','reduced_patches');

end
