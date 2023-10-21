function [centers,idxs] = anderbergy(data, numOfClusters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numOfSamples, numOfDims] = size(data);

if numOfClusters > numOfSamples
   disp('Error: Number of Cluster > Number of Samples.');
   idxs = [];
   centers = [];
   return;
end;

centers = data(1:numOfClusters,:);
idxs = [1:numOfClusters];

[dmi,C1,C2,MatrixD]=MinDistanceBetweenCenters(centers);

for i=(numOfClusters+1):numOfSamples
   multiPoint = ones(numOfClusters,1)*data(i,:);
   dist = sum((centers-multiPoint).^2, 2);
      
   [sortedDist, distIdx] = sort(dist);
   dj = sortedDist(1);
   cj = distIdx(1);

   if dj > dmi;
      centers_Cand1 = centers;
      centers_Cand1(C1,:) = data(i,:);  
      [dmi_Cand1,C1_Cand1,C2_Cand1,MatrixD_Cand1]=MinDistanceBetweenCenters(centers_Cand1);
      
      centers_Cand2 = centers;
      centers_Cand2(C2,:) = data(i,:);  
      [dmi_Cand2,C1_Cand2,C2_Cand2, MatrixD_Cand2]=MinDistanceBetweenCenters(centers_Cand2);
            
      if dmi_Cand1 > dmi_Cand2;
         centers = centers_Cand1;
         dmi = dmi_Cand1;
         C1 = C1_Cand1;
         C2 = C2_Cand1;
         MatrixD = MatrixD_Cand1;
         idxs(C1)=i; 
      else
         centers = centers_Cand2;
         dmi = dmi_Cand2;
         C1 = C1_Cand2;
         C2 = C2_Cand2;
         MatrixD = MatrixD_Cand2;
         idxs(C2)=i;
      end 
      
   else      		      
      minDistToRemainingCenters = sortedDist(2);
      
      fullMatrixD = MatrixD + MatrixD';
      DistCjToOthersCenters = [fullMatrixD(cj,1:(cj-1)) fullMatrixD(cj,(cj+1):size(fullMatrixD,2))];     
                  
      if ((min(DistCjToOthersCenters)) < minDistToRemainingCenters)
         centers(cj,:) = data(i,:);
         idxs(cj)=i;
         [dmi,C1,C2,MatrixD]=MinDistanceBetweenCenters(centers);
      end  
       
	end 
end 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dmi,C1,C2,matrixD]=MinDistanceBetweenCenters(centers);

numOfCenters = size(centers,1);

matrixD = zeros(numOfCenters);
d=[];
k=1;
for i=1:numOfCenters 
   for j=i+1:numOfCenters 
      matrixD(i,j)= (centers(i,:)-centers(j,:))*(centers(i,:)-centers(j,:))';
      d(k)= matrixD(i,j);
      k=k+1;
   end 
end 

[dmi, minIdx] = min(d);
[minRow, minCol] = find(matrixD==dmi);

C1 = minRow(1);
C2 = minCol(1);

return
