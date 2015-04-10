load gainChannel
i=1:100;
x=repmat(i,1,64);
y=ones(1,100);
for ii=2:64
    aa = ii .* ones(1,100);
    y=[y aa];
end

% gainChannel = zeros( TOTAL_USER, TOTAL_SUB, LOOP );
z=zeros(1,6400);
p=1;
for jj=1:64
    for j=1:100
        z(1,p)=gainChannel(j,jj,3);
        p=p+1;
    end
end

 stem3(x,y,z)
  view(0,0)