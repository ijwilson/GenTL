#include "binode.H"

inline int WithinCalculations(int L, int R,int dL,int dR)
{
  if (L==0||R==0) return 0;
  int tmp=(L==1)?1+dL+dR:0;
  if (R==1) tmp += 1+dL+dR;
  return tmp;
}

int BetweenCalculations(int CC[2],int L[2], int R[2],int dL[2],int dR[2])
{
  if (CC[0]==0||CC[0]==CC[1]) return 0;   // all Cases or Controls
  int tmp=0;
  if (L[0]==0) tmp += L[1]*(dL[1]+dR[0]+1);
  if (R[0]==R[1]) tmp += R[0]*(dL[1]+dR[0]+1);
  if (L[0]==L[1]) tmp += L[0]*(dL[0]+dR[1]+1);  
  if (R[0]==0) tmp += R[1]*(dL[0]+dR[1]+1);
            
  return tmp;
}


int BetweenCalculationsB(int CC[2], int L[2], int R[2],int dL[2],int dR[2])
{
  if (CC[0]==0||CC[0]==CC[1]) return 0;   // all Cases or Controls
  int tmp=0;
  // if (L[0]==0) tmp += L[1]*(dL[1]+dR[0]+1);
  if (R[0]==R[1]) tmp += R[0]*(dL[1]+dR[0]+1);
  if (L[0]==L[1]) tmp += L[0]*(dL[0]+dR[1]+1);  
  //if (R[0]==0) tmp += R[1]*(dL[0]+dR[1]+1);
            
  return tmp;
}

std::vector<double> binode::TreeDistanceStatistic(const std::vector<int> &Cases, int CC[2], int d[2])
{
    std::vector<double> PartialSums(4,0.0);
    // these partial sums are for minD(Case,Case), minD(Control,Control), minD(Case,Control), minD(Control,Case)
    if (isleaf()) {
      CC[0]=CountIntersection(Cases,labels);
      CC[1]=static_cast<int>(labels.size());
      d[0]=d[1]=0;
      // all sums are zero at a leaf
    } else {
      int LeftDistance[2],RightDistance[2];
      int LeftCC[2],RightCC[2];
      std::vector<double> LeftStat=left->TreeDistanceStatistic(Cases,LeftCC,LeftDistance);
      std::vector<double> RightStat=right->TreeDistanceStatistic(Cases,RightCC,RightDistance);
      
      if (LeftCC[0]>0&&RightCC[0]>0) d[0]=std::min(LeftDistance[0],RightDistance[0])+1;
      else if (LeftCC[0]>0) d[0] = LeftDistance[0]+1;
      else  if (RightCC[0]>0) d[0] = RightDistance[0]+1;
      else d[1]=9999;
    
      if (LeftCC[0]<LeftCC[1]&&RightCC[0]<RightCC[1]) d[1]=std::min(LeftDistance[1],RightDistance[1])+1;
      else if (LeftCC[0]<LeftCC[1]) d[1] = LeftDistance[1]+1;
      else if (RightCC[0]<RightCC[1])  d[1] = RightDistance[1]+1;
      else d[1]=9999;

      CC[0]=LeftCC[0]+RightCC[0];
      CC[1]=LeftCC[1]+RightCC[1];

      PartialSums[0] = LeftStat[0]+RightStat[0]+static_cast<double>(WithinCalculations(LeftCC[0],RightCC[0],LeftDistance[0],RightDistance[0]));
      PartialSums[1] = LeftStat[1]+RightStat[1]+static_cast<double>(WithinCalculations(LeftCC[1]-LeftCC[0],RightCC[1]-RightCC[0],LeftDistance[1],RightDistance[1]));
      PartialSums[2] = LeftStat[2]+RightStat[2]+static_cast<double>(BetweenCalculations(CC,LeftCC,RightCC,LeftDistance,RightDistance));
  
    }
    return PartialSums;
  }

std::pair<int,int> binode::RecurseLeftRightDistances(const std::vector<int> &Cases)
{
  distance=TNT::Array2D<int>(2,3);

  if (isleaf()) {
    localCC[0]=CountIntersection(Cases,labels);
    localCC[1]=static_cast<int>(labels.size());
    if (localCC[0]>0&&localCC[1]>localCC[0]) return std::pair<int,int>(0,0);
    else if (localCC[0]>0) return std::pair<int,int>(0,999);
    else if (localCC[1]>localCC[0]) return std::pair<int,int>(999,0);
    else return std::pair<int,int>(999,999);
    // all sums are zero at a leaf
  } else {
    std::pair<int,int> leftDistance = left->RecurseLeftRightDistances(Cases);
    std::pair<int,int> rightDistance = right->RecurseLeftRightDistances(Cases);

    distance[0][1] = leftDistance.first;
    distance[1][1] = leftDistance.second;
    distance[0][2] = rightDistance.first;
    distance[1][2] = rightDistance.second;
    localCC[0]=left->localCC[0]+right->localCC[0];
    localCC[1]=left->localCC[1]+right->localCC[1];

    return std::pair<int,int>(std::min(leftDistance.first,rightDistance.first)+1,
                              std::min(leftDistance.second,rightDistance.second)+1);
  }
}


std::vector<int>  binode::RecurseCherries(const std::vector<int> &Cases)
{

  if (labels.size()==2) { // a cherry, go no further
    std::vector<int> cherries(3,0); 
    cherries[2-CountIntersection(Cases,labels)]=1;
    return cherries;
  }  else if (isleaf()) {
    // don't know what else to do here.  For some sizes and combinations of 
    // cases and controls we can guarantee some cherries
    // I think the best thing is to split until cherries - as Mauro likes!
    return std::vector<int>(3,0);
  } else {
    std::vector<int> LeftCherries= left->RecurseCherries(Cases);
    std::vector<int> RightCherries = right->RecurseCherries(Cases);
    
    LeftCherries[0] += RightCherries[0];
    LeftCherries[1] += RightCherries[1]; 
    LeftCherries[2] += RightCherries[2];
    return LeftCherries;
  }
}

std::vector<double>  binode::RecurseSumHeight(const std::vector<int> &Cases,int CC[2])
{
  std::vector<double> SumHeights(2);
  if (isleaf()) {
    CC[0]=CountIntersection(Cases,labels);
    CC[1]=static_cast<int>(labels.size());
    return std::vector<double>(2,0.0);
  } else {
    int CCLeft[2],CCRight[2];
    std::vector<double> leftHeight = left->RecurseSumHeight(Cases,CCLeft);
    std::vector<double> rightHeight = right->RecurseSumHeight(Cases,CCRight);
    CC[0] = CCLeft[0] + CCRight[0];
    CC[1] = CCLeft[1] + CCRight[1];
    std::vector<double> heights(2,0.0);
    heights[0] = leftHeight[0] + CC[0] + rightHeight[0];
    heights[1] = leftHeight[1] + CC[1]-CC[0] + rightHeight[1];
    return heights;
  }
}
