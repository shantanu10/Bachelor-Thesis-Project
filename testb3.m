sam=load('Haberman.txt');
figure(1)
plot(sam(:,1),sam(:,2),'k*','MarkerSize',5);
title 'Haberman Data';
xlabel 'Attrinute 1';
ylabel 'Attrinute 2';
rng(1); % For reproducibility
sam = sam(:,1:2);
[S,D] = size(sam);

%D=D-1;
sam1=sam(:,1:D);
N = 20;
A = 0.5;  %loudness
r = 0.5;   %pulse rate
K = 3;
M = 20;
%frequency range
Qmin = 0;
Qmax = 0.02;

%iteration parameter
tol = 10^(-5);
N_iter = 0;

time=[];
%initial arrays
Q = zeros(N,1); %frequency
v = zeros(N,D*K); %velocities
%tic
%%
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
%  The initial locations of n bats
x=zeros(N,K*D);
cen=zeros(K,D);
fitt2=fitt;
% best solution
fg=inf;
pg=x(1,:);

%%
%----------start---------%
  tic
     for count=1:M
         count
       
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
    %update.move
    [lightn,index] = sort(fitt2);
    x=x(index,:);
    [fmin,index] = min(fitt2);
    best = x(1);
      %lightn
    
    for i=1:N
        if fitt2(i)<fg
            pg=x(i,:);
            fg=fitt2(i);
            fljg=clmat(i,:);
            bestx=i;
        end
    end
    %fljg
    fg
    checkline(count)=fg;  
     
  %while (fmin>tol)
  %count
  %loop over all bats/solutions
    for i=1:N,
        Q(i)=Qmin+(Qmin-Qmax)*rand;
        v(i,:)=v(i,:)+(x(i,:)-best)*Q(i);
        x(i,:)=x(i,:)+v(i,:);
        %pulse rate
        if rand>r
            x(i,:)=best+0.01*randn(1,D*K);
        % update centriod
                for z=1:K
                    cen(z,:)=x(i,(z-1)*D+1:z*D);
                end
        end
        %evaluate new solutions
        %Fnew = Fun(S(i,:));
                for p=1:S
                    tmp1=zeros(1,K);
                    for k=1:K
                        %the distance between sample and central point
                        tmp1(k)=sum((sam1(p,:)-cen(k,:)).^2);
                    end
                    %clustering
                    [tmp2 clmat(i,p)]=min(tmp1);
                end
            
                %if the solution improves or not too loud 
                    if (tmp2<fitt2(i)&(rand<A))
                        x(i,:)=x(i,:);
                        fitt2(i)=tmp2;
                    end
                    %update the current best
                    if tmp2<=fmin
                        best=x(i,:);
                        fmin=tmp2;
                    end
         
        
    end
    N_iter=N_iter+N;
    
 % time(count)=toc
%% get precision and recall

    
    
    end
%% % get centre     
      toc
%                 for z=1:K
%                     cen(z,:)=x(bestx,(z-1)*D+1:z*D);
%                 end
    
%% inter-cluster distance
interdist=0;
% for i=1:K
%      aa=find(ww(:,i)==1);
%     for j=1:K
%         bb=find(ww(:,i)==1);
%         if i~=j
%             if length(aa)~=0 & length(bb)~=0
%                 for k=1:length(aa) 
%                     for  n=1:length(bb) 
%                      interdist=interdist+sum((sam1(aa(k),:)-sam1(bb(n),:)).^2);
%                     end
%                 end
%             end
%         end
%      end
% end


for i=1:K
    for j=1:K
        if i~=j
            interdist=interdist+sum((centmp(i,:)-centmp(j,:)).^2);
        end
    end
end
interdist=interdist/2;
%
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

%%
%time=toc
figure(2)
plot(checkline);
title(' bat clustering convergence');
figure(3)
%plot3(cen(:,1),cen(:,2),cen(:,3),'o');grid;box
plot(cen(:,1),cen(:,2),'o');grid;box
title('bat clustering')
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
   legend('center','cluster 1','cluster 2','cluster 3');
rotate3d
std1=std2(sam1(a,:));
std22=std2(sam1(b,:));
std3=std2(sam1(c,:));
std4=std2(sam1(d,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quality of clusters %%%%%%%%%%%%%%%%%%%%%%
data = sam1;
ids = fljg;
%%%%%% (1) calculation of Silhoutte coefficient%%%%%
figure(9)
[silh,h] = silhouette(data, ids);
avrgScore = mean(silh);
fprintf('Mean value of the silhoutte value over all the points is %f\n',avrgScore);

%%%% (2)calculation of Dunn Index%%%%%%%%
distM=squareform(pdist(data));
disp(sprintf('Dunns index for kmeans %d', dunns(K,distM,ids)));
figure(10)
scatter(data(:,1),data(:,2),10,ids,'filled');title('Clusters found by bat algorithm');axis equal;
    