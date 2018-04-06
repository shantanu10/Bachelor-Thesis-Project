sam=load('SALES1.txt');
figure(1)
plot(sam(:,1),sam(:,2),'k*','MarkerSize',5);
title 'Fisher''s Iris Data';
xlabel 'Petal Lengths (cm)';
ylabel 'Petal Widths (cm)';
rng(1); % For reproducibility


N=20;       
% randomness 0-1
alpha=0.2;   
% Absorption coefficient
gamma=1.0;   
% number of iterations
M=20;       
% number of clusters
K=2;         
% number of samples and number of attributes
[S,D]=size(sam);
%D=D-1;
sam1=sam(:,1:D);
% initial light intensity
%lightn=rand(N,K*D);
time=[]
%plot3(sam(:,1),sam(:,2),sam(:,3),'*');grid;box;
  % line(sam(:,1),sam(:,2),sam(:,3),'linestyle','none','marker','*','color','g');
t=1;
% initialize clustering matrix clamt is a N*S
for i=1:N
    clmat(i,:)=randperm(S);
    clmat(i,clmat(i,:)>K)=ceil(rand(1,sum(clmat(i,:)>K))*K);
end

%%
% initialize best fitness
fitt=inf*ones(1,N);
% best clustring
fljg=clmat(1,:);
%  The initial locations of n fireflies
x=zeros(N,K*D);
cen=zeros(K,D);
fitt2=fitt;
% best solution
fg= inf;
pg=x(1,:);

%%
%----------start---------%
while t<M
    tic
   % t
    for i=1:N
        ww=zeros(S,K);
        % if belong class, the value is 1 otherwise 0
        for ii=1:S
            ww(ii,clmat(i,ii))=1;
        end
        ccc=[];tmp=0;
        for j= 1:K
            % sum of attributes values
            sumcs = sum(ww(:,j)*ones(1,D).*sam1);
            % number of attributes
            countcs = sum(ww(:,j));
            % class center
            % attributes' number is 0, the center is 0
            if countcs==0
                cen(j,:)=zeros(1,D);
            % otherwise the center is mean of attributes value
            else
                cen(j,:)=sumcs/countcs;
                 centmp(j,:)=cen(j,:);
            end
            ccc=[ccc,cen(j,:)];
            aa=find(ww(:,j)==1);
            if length(aa)~=0
                for k=1:length(aa)
                    tmp=tmp+(sum((sam1(aa(k),:)-cen(j,:)).^2));
                end
            end
        end
        x(i,:)=ccc;
        %best fitness
        fitt2(i)=tmp;    
    end
    % update.move
    [lightn,index] = sort(fitt2);
    x=x(index,:);
      
    for i=1:N
        if fitt2(i)<fg
            pg=x(i,:);
            fg=fitt2(i);
            fljg=clmat(i,:);
        end
    end

    fg
    %fljg
    checkline(t)=fg;
    
    for i=1:N
        for j=1:N
            r=0;
            for k=1:K*D
            r=sqrt((x(1,k)-x(i,k))^2+r);
            end

            % the firefly movement
            if lightn(i)<lightn(j)
                beta0=1;
                beta=beta0*exp(-gamma*r.^2);
                x(i,:)=x(i,:).*(1-beta)+x(j,:).*beta+alpha.*(rand-0.5);
             % update centriod
                for z=1:K
                    cen(z,:)=x(i,(z-1)*D+1:z*D);
                end
             %update
             %getting new central point
                for p=1:S
                    tmp1=zeros(1,K);
                    for k=1:K
             %the distance between sample and central point
                    tmp1(k)=sum((sam1(p,:)-cen(k,:)).^2);
                    end
             %clustering
                    [tmp2 clmat(i,p)]=min(tmp1);
               end


            end
       end
        %fljg=clmat(i,:);
    end
        %lightn
    time(t)=toc;
t=t+1;
       % checkline(t)=fg;
end
%%
interdist=0;
for i=1:K
    for j=1:K
        if i~=j
        interdist=interdist+sum((centmp(i,:)-centmp(j,:)).^2);
        end
    end
end
interdist=interdist/2;

    
%%
%display
a=[];
b=[];
c=[];
d=[];
for i=1:S
    if fljg(i)==1.0
     a(end+1)=i;
    elseif fljg(i)==2.0
        b(end+1)=i;
    elseif fljg(i)==3.0
        c(end+1)=i;
    elseif fljg(i)==4.0
        d(end+1)=i;
    end
end

     
    
    
    
    
    
    toc
    
    

    
%%       Abandoned
    
    
    
%% get best solution
    for i=1:N
         if fitt2(i)<fg
             pg=x(i,:);
             fg=fitt2(i);
             %current best clutering
             fljg=clmat(i,:);
         end
    end
fprintf('Initial centroids should be:');    
fprintf( '%d, %d \n' ,pg);   


figure(2)
plot(checkline);
title(' firefly clustering convergence');
figure(3)
%plot3(cen(:,1),cen(:,2),cen(:,3),'o');grid;box
plot(cen(:,1),cen(:,2),'ko');grid;box
title('firefly_clustering')
xlabel('First attribute')
ylabel('Second attribute')
YY=[1 2 3 4];
% index1 = find(YY(1) == best_solution)
% index2 = find(YY(2) == best_solution)
% index3 = find(YY(3) == best_solution)
% index4 = find(YY(4) == best_solution)
   line(sam(a,1),sam(a,2),'linestyle','none','marker','*','color','g');
   line(sam(b,1),sam(b,2),'linestyle','none','marker','*','color','r');
   line(sam(c,1),sam(c,2),'linestyle','none','marker','+','color','b');
   line(sam(d,1),sam(d,2),'linestyle','none','marker','s','color','b');
   legend('center','cluster 1','cluster 2','cluster 3')
rotate3d
std1=std2(sam1(a,:));
std22=std2(sam1(b,:));
std3=std2(sam1(c,:));
std4=std2(sam1(d,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quality of clusters %%%%%%%%%%%%%%%%%%%%%%
data = sam1;
ids = fljg;
%%%%%% (1) calculation of Silhoutte coefficient%%%%%
figure(4)
[silh,h] = silhouette(data, ids);
avrgScore = mean(silh);
fprintf('Mean value of the silhoutte value over all the points is %f\n',avrgScore);

%%%% (2)calculation of Dunn Index%%%%%%%%
distM=squareform(pdist(data));
disp(sprintf('Dunns index is %d', dunns(K,distM,ids)));
figure(5)
scatter(data(:,1),data(:,2),10,ids,'filled');title('Clusters found by firefly algorithm');axis equal;

