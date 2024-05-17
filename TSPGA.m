clear all
close all hidden
x = [1 2 5 7 8 10 12 13 14 15];
y = [4 14 2 10 6 12 9 3 6 18];

A = zeros(size(x,2),size(y,2));
for i = 1:size(x,2)
    for j = 1:size(y,2)
        A(i,j) = norm([x(i),y(i)]-[x(j),y(j)]);
    end
end

% for j = 1:10
%     cntSuc = 0;
%     for i = 1:250
        [tbst, lenbst, fst, bestAvgWorstMat, itrCnt, stat] = runTSPGA(A,100,80,0.6);
%         if getDist(A, tbst) < 60
%             cntSuc = cntSuc + 1;
%         end
%     end
%     perSuc(j)=cntSuc/250;
% end
% perSuc


figure(1)
hold on
labels = {'A  ','B  ','C  ','D  ','E  ','F  ','G  ','H  ','I  ','J  '};
plot(x,y,"o",'Linewidth',2)
text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

plot(x(tbst), y(tbst),"c",'Linewidth',1.5)
plot([x(tbst(1)) x(tbst(end))], [y(tbst(1)) y(tbst(end))], "c",'Linewidth',1.5)

plot(x(fst), y(fst),'k--')
plot([x(fst(1)) x(fst(end))], [y(fst(1)) y(fst(end))],'k--')

xlim([0 16]); ylim([0 19]);
xlabel("x"); ylabel("y");
hold off

figure(2)
hold on
plot(1:itrCnt,bestAvgWorstMat(1,1:itrCnt),'r','Linewidth',1)
plot(1:itrCnt,bestAvgWorstMat(2,1:itrCnt),'b','Linewidth',1)
plot(1:itrCnt,bestAvgWorstMat(3,1:itrCnt),'k','Linewidth',1)
xlabel("Iteration"); ylabel("Tour Length");
legend(["best", "average", "worst"]);
hold off


function [tbst, lenbst, fst, bestAvgWorstMat, itrCnt, stat] = runTSPGA(A, itrBnd,pSz,p)
%TSPGA solves a given traveling sales problem with a genetic alogrithm 
%   [tbst, lenbst ,itrCnt , stat] = TSPGA(A, itrBnd,pSz ,p) takes in A, an
%   n x n matiex such that A_ij is the distance between cities i and j.
%   itrBound an optional iteration boumd with default value 10n. psZ an
%   optional parameter that gives the number of tours in a population, the
%   defult value is 4n. p an optional parameter that gives the maxumum
%   mutation percentage, the defult value is 0.4. TSPGA returns four
%   parameter tbst the best tour of the entire search, which may not be in
%   the final population due to mutation. lenbst the best tour of the
%   entire search. itrCnt the number of iteration performed. stat a status
%   variable.
%
%   stat = 0: The last population contains the best tour.
%   stat = 1: The last population does not contain the best tour.
arguments
    A
    itrBnd = 10*size(A,1)
    pSz = 4*size(A,1)
    p = 0.4
end


% A = [0 10.1 4.5 8.5 7.3; 10.1 0 12.4 6.4 10; 4.5 12.4 0 8.2 5; 8.5 6.4 8.2 0 4.1; 7.3 10 5 4.1 0;]

numCities = size(A,1);
tours = zeros(numCities, pSz);

s = RandStream('mlfg6331_64'); 
for i=1:pSz
    tours(:,i) = (datasample(s,1:numCities,numCities,'Replace',false))';
end

%intial variables
tbst = tours(:,1);
fst = tbst;
lenbst = tbst;
itrCnt = 0;
stat = 1;
bestDist = Inf;
noImpCount = 0;

bestAvgWorstMat = zeros(3,itrBnd);

% tours
for i=1:itrBnd
    itrCnt = i;

    %get fitness
    dists = zeros(1,pSz);
    for j=1:pSz
        tempDist = 0;
        for k=1:numCities-1
            tempDist = tempDist + A(tours(k,j),tours(k+1,j));
        end
        tempDist = tempDist + A(tours(1,j),tours(end,j));
        dists(j) = tempDist;
    end

    %reorder based on fitness, most fit 1   at col 1
    [dists, distsOrder] = sort(dists);
    tours = tours(:,distsOrder);

    bestAvgWorstMat(:,itrCnt) = [dists(1) mean(dists) dists(end)]';

    % test if curr top tour is better then the best
    if dists(1)<bestDist
        bestDist = dists(1); 
        lenbst = tours(:,1);
    end

    tbst = tours(:,1);    

    %crossover
    tempTours = zeros(numCities, pSz);
    crossSize = ceil(pSz/2);
    for j=1:pSz
        % best tour mates with itself
        if j==1
            par1 = tours(:,1);
            par2 = tours(:,1);
            par1val = randi(numCities,1);
            par1Pos = find(par1==par1val,1);
            par2([par1Pos,find(par2==par1val,1)]) = par2([find(par2==par1val,1),par1Pos]);
            tempTours(:,j) = par2;
        % best 50% mate randomly with each other
        else
            par1 = tours(:,randi(crossSize,1));
            par2 = tours(:,randi(crossSize,1));
            par1val = randi(numCities,1);
            par1Pos = find(par1==par1val,1);
            par2([par1Pos,find(par2==par1val,1)]) = par2([find(par2==par1val,1),par1Pos]);
            tempTours(:,j) = par2;
        end
    end
    tours = tempTours;


    %mutate
    for c=1:pSz
        if rand(1)<p
            curTour = tours(:,c);
            r1 = randi(numCities,1);
            r2 = randi(numCities,1);
            curTour([r1,r2]) = curTour([r2,r1]);
            tours(:,c) = curTour;
        end
    end

    if tbst == lenbst
        noImpCount = noImpCount + 1;
    else
        noImpCount = 0;
    end

    % return if no change in 2n iterations
    if noImpCount >= 2*size(A,1), return; end

end

finalTopDist = 0;
for k=1:numCities-1
    finalTopDist = finalTopDist + A(tours(k,1),tours(k+1,1));
end

if finalTopDist==bestDist, stat = 0; end

end

function dist = getDist(A, tour)
    dist = 0;
    for k=1:size(A,1)-1
        dist = dist + A(tour(k),tour(k+1));
    end
    dist = dist + A(tour(1),tour(end));
end
