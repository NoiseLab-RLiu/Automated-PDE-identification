function [CVERR] = CrossValid5(Xs,y)
%5-fold cross validation, CVERR is the validation error, with first row
%corresponding to sparsity=0.
y = y/norm(y);
dict = Xs;
length = size(y,1)/5;

randseed = 1;
rng(randseed)
r = randperm(size(y,1));%randperm
length = size(r,2)/5; %% PAY SPECIAL NOTICE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dictfold1=dict(r(0*length+1:1*length),:); %Expriments suggests that do not use random selection here. 
dictfold2=dict(r(1*length+1:2*length),:);
dictfold3=dict(r(2*length+1:3*length),:);
dictfold4=dict(r(3*length+1:4*length),:);
dictfold5=dict(r(4*length+1:5*length),:);

yfold1=y(r(0*length+1:1*length));
yfold2=y(r(1*length+1:2*length));
yfold3=y(r(2*length+1:3*length));
yfold4=y(r(3*length+1:4*length));
yfold5=y(r(4*length+1:5*length));


dicv1=[dictfold2;dictfold3;dictfold4;dictfold5];
dicv2=[dictfold1;dictfold3;dictfold4;dictfold5];
dicv3=[dictfold1;dictfold2;dictfold4;dictfold5];
dicv4=[dictfold1;dictfold2;dictfold3;dictfold5];
dicv5=[dictfold1;dictfold2;dictfold3;dictfold4];

ycv1 = [yfold2;yfold3;yfold4;yfold5];
ycv2 = [yfold1;yfold3;yfold4;yfold5];
ycv3 = [yfold1;yfold2;yfold4;yfold5];
ycv4 = [yfold1;yfold2;yfold3;yfold5];
ycv5 = [yfold1;yfold2;yfold3;yfold4];

dictfold = cell(5,1);
dictfold{1} = dictfold1;
dictfold{2} = dictfold2;
dictfold{3} = dictfold3;
dictfold{4} = dictfold4;
dictfold{5} = dictfold5;

yfold = cell(5,1);
yfold{1} = yfold1;
yfold{2} = yfold2;
yfold{3} = yfold3;
yfold{4} = yfold4;
yfold{5} = yfold5;

dicv = cell(5,1);
dicv{1} = dicv1;
dicv{2} = dicv2;
dicv{3} = dicv3;
dicv{4} = dicv4;
dicv{5} = dicv5;

ycv = cell(5,1);
ycv{1} = ycv1;
ycv{2} = ycv2;
ycv{3} = ycv3;
ycv{4} = ycv4;
ycv{5} = ycv5;

CVERR=zeros(size(Xs,2)+1,6);
for fold=1:5
error=[];
for i=1:size(Xs,2)
    [x,error_abd] = OMP_N(dicv{fold},ycv{fold},i); % train
    error=[error,norm(yfold{fold}-dictfold{fold}*x)^2] % test error
end

error = [norm(yfold{fold})^2;error']; % 1st entry is for sparsity 0.

CVERR(:,fold) = error;
end

CVERR(:,6)=5*mean(CVERR(:,1:5),2);
end

