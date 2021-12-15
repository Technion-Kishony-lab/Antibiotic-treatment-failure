function c = my_fit(X,Y,n0,link)

if nargin<3
    n0 = 0 ; % minimal number of zeros and ones
end
if nargin<4
    link = 'logit' ;
end
sumX0 = sum(X<max(X)) ;
sumX1 = sum(X>min(X)) ;
kv = sumX0>=n0 & sumX1>=n0 ; % variables to include
X = X(:,kv) ;

[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial',link) ;
c.coef = glm.stats.beta ;
c.se = glm.stats.se ;
c.p = glm.stats.p ;
c.t = glm.stats.t ;
c.dfe = glm.stats.dfe ;
c.covb = glm.stats.covb ;
c.rmse = glm.stats.sfit ;
c.coeffcorr = glm.stats.coeffcorr ;

kv = [true, kv] ;
c.coef = addnans(c.coef,kv,1) ;
c.se   = addnans(c.se  ,kv,1) ;
c.p    = addnans(c.p   ,kv,1) ;
c.t    = addnans(c.t   ,kv,1) ;
c.covb = addnans(c.covb,kv,3) ;
c.coeffcorr = addnans(c.coeffcorr,kv,3) ;

c.sumX0 = sumX0 ;
c.sumX1 = sumX1 ;
c.sumY0 = sum(Y==0) ;
c.sumY1 = sum(Y>0) ;
c.totallines = size(X,1) ;
c.remove = ~kv ;

end

function xn = addnans(x,k,dim)
switch dim
    case 1
        xn(k,:) = x ;
        xn(~k,:) = nan ;
    case 2
        xn(:,k) = x ;
        xn(:,~k) = nan ;
    case 3
        xn(k,k) = x ;
        xn(:,~k) = nan ;
        xn(~k,:) = nan ;
end
end

