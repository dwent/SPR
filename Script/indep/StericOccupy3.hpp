//Version 2
//Date 2008.2.29.
//Author caoyang
// g++ StericOccupy.cpp -o StericOccupy
//Designed for search neighbour
//Date 2008.7.8
//latticesize  设置太小的话会使得格点数目超过整型数据最大范围或者内存不够用而报错
//增加VertexAmount(vector<int> &effpoint, vector<vector<int> > &index);
//区分单个原子的值
//2009.3.5.
#ifndef Steric3Head
#define Steric3Head

#include <vector>
#include <cstdlib>

using namespace std;


//构造坐标vector<float> &xyz，得到半径radii的坐标点
int ShellBuilding(vector<float> &xyz, float radii, float step, vector<vector<float> > &group)
{
      int i, j, k, num=0, max=int(radii/step)+1, min=max/2;
      float up=(radii/step+0.5)*(radii/step+0.5), down=(radii/step-0.5)*(radii/step-0.5);
      //cout<<"Max "<<max<<" Min "<<min<<endl;
      vector<float> tmp(3, 0);
      for(i=-max; i<=max; ++i)
      {
         for(j=-max; j<=max; ++j)
         {
            for(k=-max; k<=max; ++k)
            {
               if(i*i+j*j+k*k> up|| i*i+j*j+k*k< down) continue;
               tmp[0]=xyz[0]+step*i;
               tmp[1]=xyz[1]+step*j;
               tmp[2]=xyz[2]+step*k;
               
               group.push_back(tmp);
               num++;
               //cout<<"ATM  "<<tmp[0]<<"  "<<tmp[1]<<"  "<<tmp[2]<<endl;
            }
         }
      }
      return num;
}

//构造坐标vector<float> &xyz，得到半径radii内的的坐标点
int SphereBuilding(vector<float> &xyz, float radii, float step, vector<vector<float> > &group)
{
      int i, j, k, num=0, max=int(radii/step)+1, min=max/2;
      float up=(radii*radii)/(step*step);
      //cout<<"Max "<<max<<" Min "<<min<<endl;
      vector<float> tmp(3, 0);
      for(i=-max; i<=max; ++i)
      {
         for(j=-max; j<=max; ++j)
         {
            for(k=-max; k<=max; ++k)
            {
               if(i*i+j*j+k*k> up) continue;
               tmp[0]=xyz[0]+step*i;
               tmp[1]=xyz[1]+step*j;
               tmp[2]=xyz[2]+step*k;
               
               group.push_back(tmp);
               num++;
               //cout<<"ATM  "<<tmp[0]<<"  "<<tmp[1]<<"  "<<tmp[2]<<endl;
            }
         }
      }
      return num;
}

//构造坐标vector<float> &xyz，得到半径radii内的的坐标点,及其势能值potential
int PotentialSphereBuilding(vector<float> &xyz, float radii, float step, float potential_rate, vector<vector<float> > &group, vector<float> &potential)
{
      int i, j, k, num=0, max=int(radii/step)+1, min=max/2;
      float up=(radii*radii)/(step*step);
      //cout<<"Max "<<max<<" Min "<<min<<endl;
      vector<float> tmp(3, 0);
      for(i=-max; i<=max; ++i)
      {
         for(j=-max; j<=max; ++j)
         {
            for(k=-max; k<=max; ++k)
            {
               if(i*i+j*j+k*k> up) continue;
               tmp[0]=xyz[0]+step*i;
               tmp[1]=xyz[1]+step*j;
               tmp[2]=xyz[2]+step*k;
               //if(i==0&&j==0&k==0)potential.push_back(potential_rate);
               //else potential.push_back(potential_rate/sqrt((float)(i*i+j*j+k*k)));
               potential.push_back(potential_rate*(1-sqrt((float)(i*i+j*j+k*k))/max));
               group.push_back(tmp);
               num++;
               //cout<<"ATM  "<<tmp[0]<<"  "<<tmp[1]<<"  "<<tmp[2]<<endl;
            }
         }
      }
      return num;
}
//Author caoyang
class StericOccupy
{
      public:
      StericOccupy(float latticesize, vector<vector<float> > &group1)
      {GenerateLattice(latticesize, group1);};
      StericOccupy(float latticesize, vector<int> &max3, vector<int> &seqn);
      
      //Generate Lattice 建立格子模型的格点被占据的数量表，对应原子索引
      int GenerateLattice(float latticesize, vector<vector<float> > &group1);
      //OccupyVertex 查询xyz坐标处被占据的数量
      int OccupyVertex(vector<float> &xyz);
      //NeighborIndex 查询xyz坐标附近（格子深度link）被占据的数量
      int NeighborIndex(vector<float> &xyz, int link);
      int NeighborIndex(int seqn, int link);
      //NeighborIndex 查询xyz坐标附近（半径为radii）被占据的数量
      int NeighborIndex(vector<float> &xyz, float radii);
      //FreeSpace 查询xyz坐标附近（半径为radii）空出的空间点数的数量（包括边界点外）
      int FreeSpace(vector<float> &xyz, float radii);
      //FreeSpace 查询xyz坐标附近（半径为radii）空出的空间点数（包括边界点外）占有该空间的比例
      float FreeSpaceRatio(vector<float> &xyz, float radii);
      float FreeSpaceRatio(int seqn, float radii);
      //OccupyIndex 查询xyz坐标处被占据的数量, 对应原子索引index
      int OccupyIndex(vector<float> &xyz, vector<int> &index);
      //NeighborIndex 查询xyz坐标附近（格子深度link）被占据的数量, 对应原子索引NeiIndex
      int NeighborIndex(vector<float> &xyz, int link, vector<int> &NeiIndex);
      //BoundaryFaceIndex识别边界面的点
      int BoundaryFaceIndex(vector<float> &xyz, int link);
      int BoundaryFaceIndex(vector<float> &xyz);
      int BoundaryFaceIndex(int seqn);
      //GeneratePotential添加格点的势能
      int GeneratePotential(vector<vector<float> > &group1, vector<float> &potential);
      //VertexPotential 查询xyz坐标处的势能
      float VertexPotential(vector<float> &xyz);
      float VertexPotential(int seqn);
      int ShowSize();
      //VertexAmount 查询所有占据点的总数量，对应了空间的体积量
      int VertexAmount(vector<int> &effpoint);
      int VertexAmount(vector<int> &effpoint, vector<vector<int> > &index);
      //测试所用
      void IdExchange(int index);
      void CopyLattice(vector<int> &max3){max3=max;};
      
      private:
      vector<vector<int> > lattice_index; //各个格点被占据的对应原子索引
      vector<int> lattice_vertex;         //格点被占据的数量表，整个格子模型的核心！
      vector<int> max;      //三个维度的最大值，即各个维度占据的格点数
      vector<float> center; //group1中心坐标
      vector<float> shift;  //lattice_vertex[0]的坐标
      vector<float> p_vertex; //格点的势能
      float lattice_size;   //单位立方格子边长
      int group1_num;       //整个格点模型中的原子数目

};

StericOccupy::StericOccupy(float latticesize, vector<int> &max3, vector<int> &seqn)
{
      lattice_vertex.assign(max3[0]*max3[1]*max3[2], 0);
      lattice_size=latticesize;
      if(seqn.size()>=max3[0]*max3[1]*max3[2])
      {
         cout<<"Error in Parameter Size\n";
         exit(0);
      }
      for(int i=0; i<seqn.size(); i++)
      {
         lattice_vertex[seqn[i]]=1;
      }
      max=max3;
}
//Author caoyang
int StericOccupy::GenerateLattice(float latticesize, vector<vector<float> > &group1)
{
      lattice_size=latticesize;
      group1_num=group1.size();
      int i=0, j=0;
      string atom;
      vector<int> tmp(3);
      vector<float> xyz1, xyz2;
      center.assign(3, 0);
      max.assign(3, 0);
      
      xyz1=group1[0]; xyz2=group1[0];
      for(i=0; i<group1_num; ++i)
      {
         for(j=0; j<3; ++j)
         {
            if(xyz1[j]>group1[i][j]) xyz1[j]=group1[i][j];
            if(xyz2[j]<group1[i][j]) xyz2[j]=group1[i][j];
            center[j]+=group1[i][j];
         }
      }
      shift=xyz1;
      for(j=0; j<3; ++j) 
      {
         xyz2[j]-=xyz1[j];
         max[j]=int(0.5 + xyz2[j]/lattice_size)+1;
         center[j]/=group1_num;
      }
      //cout<<"Vertex "<<max[0]*max[1]*max[2]<<endl;
      lattice_vertex.assign(max[0]*max[1]*max[2], 0);
      lattice_index.resize(max[0]*max[1]*max[2]);
      
      for(i=0; i<group1_num; ++i)
      {
         for(j=0; j<3; ++j) 
         {
            tmp[j]=int(0.5 + (group1[i][j]-shift[j])/lattice_size);
         }
         lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]]+=1;
         lattice_index[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]].push_back(i);
      }
      //for(i=0; i<lattice_vertex.size(); ++i)
      {
         //cout<<lattice_index[i].size()<<"  ";
      }
      //cout<<endl;
      return lattice_vertex.size();
}
int StericOccupy::OccupyVertex(vector<float> &xyz)
{
      vector<int> tmp(3);
      for(int j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
         if(tmp[j]>=max[j] || tmp[j]<0) 
         {
            cout<<tmp[j]<<" <=> "<<max[j]<<" Out of Box\n";
            return 0;
         }
      }
      return lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]];
}
/*
int StericOccupy::NeighborIndex(vector<float> &xyz, int link)
{
      int i, j, k, num=0;
      vector<int> occupyindex, tmp(3);
      //OccupyIndex(xyz, occupyindex);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
      }
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
               num+=lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k];
               //cout<<atom<<"  "<<tmp[0]+i<<"  "<<tmp[1]+j<<"  "<<tmp[2]+k<<endl;
            }
         }
      }
      return num;
}
*/
int StericOccupy::NeighborIndex(vector<float> &xyz, int link)
{
      int num=0;
      register int i, j, k;
      register int xp=int(0.5 + (xyz[0]-shift[0])/lattice_size);
      register int yp=int(0.5 + (xyz[1]-shift[1])/lattice_size);
      register int zp=int(0.5 + (xyz[2]-shift[2])/lattice_size);

      for(i=-link; i<=link; ++i)
      {
         if(xp+i<0 || xp+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(yp+j<0 || yp+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(zp+k<0 || zp+k>=max[2]) continue;
               num+=lattice_vertex[(xp+i)*max[1]*max[2]+(yp+j)*max[2]+zp+k];
            }
         }
      }
      return num;
}

int StericOccupy::NeighborIndex(int seqn, int link)
{
      int num=0;
      register int i, j, k;
      register int xp=int(seqn/(max[1]*max[2]));
      register int yp=int(seqn%(max[1]*max[2])/max[2]);
      register int zp=(seqn%(max[1]*max[2]))%max[2];
      //IdExchange(seqn);

      for(i=-link; i<=link; ++i)
      {
         if(xp+i<0 || xp+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(yp+j<0 || yp+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(zp+k<0 || zp+k>=max[2]) continue;
               num+=lattice_vertex[(xp+i)*max[1]*max[2]+(yp+j)*max[2]+zp+k];
            }
         }
      }
      return num;
}

/*
int StericOccupy::NeighborIndex(vector<float> &xyz, float radii)
{
      int i, j, k, num=0, link=int(radii/lattice_size+0.5);
      vector<int> occupyindex, tmp(3);
      //OccupyIndex(xyz, occupyindex);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
      }
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
               if(i*i+j*j+k*k>link*link) continue;
               num+=lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k];
               //cout<<atom<<"  "<<tmp[0]+i<<"  "<<tmp[1]+j<<"  "<<tmp[2]+k<<endl;
            }
         }
      }
      return num;
}
*/
//maybe there is some problem in radii. please look at FreeSpaceRatio
int StericOccupy::NeighborIndex(vector<float> &xyz, float radii)
{
      int num=0, link=int(radii/lattice_size+0.5);
      int link2=link*link;
      register int i, j, k;
      register int xp=int(0.5 + (xyz[0]-shift[0])/lattice_size);
      register int yp=int(0.5 + (xyz[1]-shift[1])/lattice_size);
      register int zp=int(0.5 + (xyz[2]-shift[2])/lattice_size);
      for(i=-link; i<=link; ++i)
      {
         if(xp+i<0 || xp+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(yp+j<0 || yp+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(zp+k<0 || zp+k>=max[2]) continue;
               if(i*i+j*j+k*k>link2) continue;
               num+=lattice_vertex[(xp+i)*max[1]*max[2]+(yp+j)*max[2]+zp+k];
               //cout<<atom<<"  "<<tmp[0]+i<<"  "<<tmp[1]+j<<"  "<<tmp[2]+k<<endl;
            }
         }
      }
      return num;
}


//与NeighborIndex(vector<float> &xyz, float radii)不同主要是统计不占空的空间点数目
int StericOccupy::FreeSpace(vector<float> &xyz, float radii)
{
      int i, j, k, num=0, link=int(radii/lattice_size+0.5);
      vector<int> tmp(3);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
      }
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) {num+=(2*link+1)*(2*link+1); continue;}
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) {num+=(2*link+1); continue;}
            for(k=-link; k<=link; ++k)
            {
               if(tmp[2]+k<0 || tmp[2]+k>=max[2]) {num++; continue;}
               if(i*i+j*j+k*k>link*link) continue;
               if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==0)
               num++;
            }
         }
      }
      return num;
}

/*
float StericOccupy::FreeSpaceRatio(vector<float> &xyz, float radii)
{
      int i, j, k, num=0, full=0, link=int(radii/lattice_size+0.5);
      vector<int> tmp(3);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
      }
      for(i=-link; i<=link; ++i)
      {
         //if(tmp[0]+i<0 || tmp[0]+i>=max[0]) {num+=(2*link+1)*(2*link+1); continue;}
         for(j=-link; j<=link; ++j)
         {
            //if(tmp[1]+j<0 || tmp[1]+j>=max[1]) {num+=(2*link+1); continue;}
            for(k=-link; k<=link; ++k)
            {
               //if(tmp[2]+k<0 || tmp[2]+k>=max[2]) {num++; continue;}
               if(i*i+j*j+k*k>link*link) continue;
               if(tmp[0]+i<0 || tmp[0]+i>=max[0] || tmp[1]+j<0 || tmp[1]+j>=max[1]
               || tmp[2]+k<0 || tmp[2]+k>=max[2])
               {
                  num++;
               }
               else if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==0)
               {
                  num++;
                  
               }
               else full++;
            }
         }
      }
      return (float)num/(num+full);
}
*/
//优化上述代码
float StericOccupy::FreeSpaceRatio(vector<float> &xyz, float radii)
{
      int num=0, full=0, link=int(radii/lattice_size+0.5);
      int link2=link*link;
      register int i, j, k;
      register int xp=int(0.5 + (xyz[0]-shift[0])/lattice_size);
      register int yp=int(0.5 + (xyz[1]-shift[1])/lattice_size);
      register int zp=int(0.5 + (xyz[2]-shift[2])/lattice_size);
      for(i=-link; i<=link; ++i)
      {
         for(j=-link; j<=link; ++j)
         {
            for(k=-link; k<=link; ++k)
            {
               if(i*i+j*j+k*k>link2) continue;
               if(xp+i<0 || xp+i>=max[0] || yp+j<0 || yp+j>=max[1]
               || zp+k<0 || zp+k>=max[2])
               {
                  num++;
               }
               else if(lattice_vertex[(xp+i)*max[1]*max[2]+(yp+j)*max[2]+zp+k]==0)
               {
                  num++;
               }
               else full++;
               
            }
         }
      }
      return (float)num/(num+full);
}
void StericOccupy::IdExchange(int index)
{
      int xp=int(index/(max[1]*max[2]));
      int yp=int(index%(max[1]*max[2])/max[2]);
      int zp=(index%(max[1]*max[2]))%max[2];
      
      int id=xp*max[1]*max[2]+yp*max[2]+zp;
      if(id!=index)
      cout<<"Error! IdExchange Test\n";
}

float StericOccupy::FreeSpaceRatio(int seqn, float radii)
{
      int num=0, full=0, link=int(radii/lattice_size+0.5);
      int link2=link*link;
      register int i, j, k;
      register int xp=int(seqn/(max[1]*max[2]));
      register int yp=int(seqn%(max[1]*max[2])/max[2]);
      register int zp=(seqn%(max[1]*max[2]))%max[2];
      //IdExchange(seqn);
      for(i=-link; i<=link; ++i)
      {
         for(j=-link; j<=link; ++j)
         {
            for(k=-link; k<=link; ++k)
            {
               if(i*i+j*j+k*k>link2) continue;
               if(xp+i<0 || xp+i>=max[0] || yp+j<0 || yp+j>=max[1]
               || zp+k<0 || zp+k>=max[2])
               {
                  num++;
               }
               else if(lattice_vertex[(xp+i)*max[1]*max[2]+(yp+j)*max[2]+zp+k]==0)
               {
                  num++;
               }
               else full++;
               
            }
         }
      }
      return (float)num/(num+full);
}


int StericOccupy::OccupyIndex(vector<float> &xyz, vector<int> &index)
{
      vector<int> tmp(3);
      for(int j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
         if(tmp[j]>=max[j] || tmp[j]<0) 
         {
            cout<<tmp[j]<<" <=> "<<max[j]<<" Out of Box\n";
            return 0;
         }
      }
      index=lattice_index[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]];
      return lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]];
}

int StericOccupy::NeighborIndex(vector<float> &xyz, int link, vector<int> &NeiIndex)
{
      int i, j, k;
      vector<int> occupyindex, tmp(3);
      //OccupyIndex(xyz, occupyindex);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
      }
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
            for(k=-link; k<=link; ++k)
            {
               if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
               NeiIndex.insert(NeiIndex.end(),lattice_index[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k].begin(), lattice_index[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k].end());
               //cout<<atom<<"  "<<tmp[0]+i<<"  "<<tmp[1]+j<<"  "<<tmp[2]+k<<endl;
            }
         }
      }
      return NeiIndex.size();
}

int StericOccupy::ShowSize()
{
   cout<<"Lattice_Size "<<max[0]<<" × "<<max[1]<<" × "<<max[2]
   <<" = "<<lattice_vertex.size()<<endl;
   return lattice_vertex.size();
}

int StericOccupy::BoundaryFaceIndex(vector<float> &xyz, int link)
{
      //不只是识别&xyz所在位置是不是边界点（界面lattice_vertex为1），怕是一些内部
      //空穴，所以还要看周围的点是不是也在边界上
      int i, j, k, num=0;
      vector<int> tmp(3);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
         if(tmp[j]>=max[j]) 
         {
            cout<<tmp[j]<<" >= "<<max[j]<<" Out of Box\n";
            return 0;
         }
      }
      k=-link;
      if(tmp[2]+k>=0 || tmp[2]+k<max[2])
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      k=link;
      if(tmp[2]+k>=0 || tmp[2]+k<max[2])
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(j=-link; j<=link; ++j)
         {
            if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      
      j=-link;
      if(tmp[1]+j>=0 || tmp[1]+j<max[1])
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(k=-link; k<=link; ++k)
         {
            if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      
      j=link;
      if(tmp[1]+j>=0 || tmp[1]+j<max[1])
      for(i=-link; i<=link; ++i)
      {
         if(tmp[0]+i<0 || tmp[0]+i>=max[0]) continue;
         for(k=-link; k<=link; ++k)
         {
            if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      
      i=-link;
      if(tmp[0]+j>=0 || tmp[0]+j<max[0])
      for(j=-link; j<=link; ++j)
      {
         if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
         for(k=-link; k<=link; ++k)
         {
            if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      i=link;
      if(tmp[0]+j>=0 || tmp[0]+j<max[0])
      for(j=-link; j<=link; ++j)
      {
         if(tmp[1]+j<0 || tmp[1]+j>=max[1]) continue;
         for(k=-link; k<=link; ++k)
         {
            if(tmp[2]+k<0 || tmp[2]+k>=max[2]) continue;
            if(lattice_vertex[(tmp[0]+i)*max[1]*max[2]+(tmp[1]+j)*max[2]+tmp[2]+k]==1)
               num++;
         }
      }
      
      if(num>2*link) return 1;
      else return 0;
   
}

int StericOccupy::BoundaryFaceIndex(vector<float> &xyz)
{
      //不只是识别&xyz所在位置是不是边界点（界面lattice_vertex为1），怕是一些内部
      //空穴，所以还要看周围的点是不是也在边界上
      int i, j, k, num=0, link=1;;
      vector<int> tmp(3);
      for(j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
         if(tmp[j]>=max[j] || tmp[j]<0) 
         {
            cout<<tmp[j]<<" >= "<<max[j]<<" Out of Box\n";
            return 0;
         }
      }
      
      if(lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]]==0) return 0;
      
      if(tmp[0]+1>=max[0]) {}
      else if(lattice_vertex[(tmp[0]+1)*max[1]*max[2]+tmp[1]*max[2]+tmp[2]]>=1) num++;
      
      if(tmp[0]-1<0) {}
      else if(lattice_vertex[(tmp[0]-1)*max[1]*max[2]+tmp[1]*max[2]+tmp[2]]>=1) num++;
      
      if(tmp[1]+1>=max[1]) {}
      else if(lattice_vertex[tmp[0]*max[1]*max[2]+(tmp[1]+1)*max[2]+tmp[2]]>=1) num++;
      
      if(tmp[1]-1<0) {}
      else if(lattice_vertex[tmp[0]*max[1]*max[2]+(tmp[1]-1)*max[2]+tmp[2]]>=1) num++;
      
      if(tmp[2]+1>=max[2]) {}
      else if(lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]+1]>=1) num++;
      
      if(tmp[2]-1<0) {}
      else if(lattice_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]-1]>=1) num++;
      
      if(num<6) return 1;
      else return 0;
}

int StericOccupy::BoundaryFaceIndex(int seqn)
{
      //不只是识别&xyz所在位置是不是边界点（界面lattice_vertex为1），怕是一些内部
      //空穴，所以还要看周围的点是不是也在边界上
      int num=0;
      vector<int> tmp(3);
      tmp[0]=int(seqn/(max[1]*max[2]));
      tmp[1]=int(seqn%(max[1]*max[2])/max[2]);
      tmp[2]=(seqn%(max[1]*max[2]))%max[2];
      //IdExchange(seqn);
      
      if(lattice_vertex[seqn]==0) return 0;
      
      if(tmp[0]+1>=max[0]) {}
      else if(lattice_vertex[seqn+max[1]*max[2]]>=1) num++;
      
      if(tmp[0]-1<0) {}
      else if(lattice_vertex[seqn-max[1]*max[2]]>=1) num++;
      
      if(tmp[1]+1>=max[1]) {}
      else if(lattice_vertex[seqn+max[2]]>=1) num++;
      
      if(tmp[1]-1<0) {}
      else if(lattice_vertex[seqn-max[2]]>=1) num++;
      
      if(tmp[2]+1>=max[2]) {}
      else if(lattice_vertex[seqn+1]>=1) num++;
      
      if(tmp[2]-1<0) {}
      else if(lattice_vertex[seqn-1]>=1) num++;
      
      if(num<6) return 1;
      
      else return 0;
}


//给格点赋予势能
int StericOccupy::GeneratePotential(vector<vector<float> > &group1, vector<float> &potential)
{
      int i=0, j=0, num=0;
      vector<int> tmp(3);
      p_vertex.assign(max[0]*max[1]*max[2], 0);
      for(i=0; i<group1.size(); ++i)
      {
         for(j=0; j<3; ++j) 
         {
            tmp[j]=int(0.5 + (group1[i][j]-shift[j])/lattice_size);
            if(tmp[j]>=max[j] || tmp[j]<0) break;
         }
         if(j==3)
         {
            p_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]]+=potential[i];
            num++;
         }
      }
      return num;
}

float StericOccupy::VertexPotential(vector<float> &xyz)
{
      vector<int> tmp(3);
      for(int j=0; j<3; ++j) 
      {
         tmp[j]=int(0.5 + (xyz[j]-shift[j])/lattice_size);
         if(tmp[j]>=max[j] || tmp[j]<0) 
         {
            cout<<tmp[j]<<" ? "<<max[j]<<"VertexPotential Out of Box\n";
            return 0;
         }
      }
      return p_vertex[tmp[0]*max[1]*max[2]+tmp[1]*max[2]+tmp[2]];
}

float StericOccupy::VertexPotential(int seqn)
{
      return p_vertex[seqn];
}


int StericOccupy::VertexAmount(vector<int> &effpoint)
{
      int i=0, amount=0;
      for(i=0; i<lattice_vertex.size(); i++)
      {
         if(lattice_vertex[i]!=0)
         {
            effpoint.push_back(i);
            amount++;
         }
      }
      return amount;
}

int StericOccupy::VertexAmount(vector<int> &effpoint, vector<vector<int> > &index)
{
      int i=0, amount=0;
      for(i=0; i<lattice_vertex.size(); i++)
      {
         if(lattice_vertex[i]!=0)
         {
            effpoint.push_back(i);
            index.push_back(lattice_index[i]);
            amount++;
         }
      }
      return amount;
}




#endif






