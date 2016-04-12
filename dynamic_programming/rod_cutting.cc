// Rod Cutting Problem
// Description:
//   Given an array of values and a rod length, find the maximum value we can get by doing the rod cutting
// Algorithm:
//   let di denotes the maximum value we can get from cutting a rod of length i, vi denotes the value of a cut of length i
//   then
//       di = 0, when i = 0;
//       di = max(di-k + vk) for 1<=k<=i, when i>=1
// Complexity:
//    let m be the number of different kinds of cut, let n be the rod length
//    time complexity: O(mn)
//    Space complexity: O(n)

#include <iostream>
#include <cmath>

using namespace std;

// output the maximum values we can get
int rod_cut(int v[], int vsize, int len){
  int values[len+1];
  values[0] = 0;
  int cur_max;
  for (int i=1; i<=len; i++){
    cur_max = -1;
    for (int j=1; j<=i; j++){
      if (j<vsize && values[i-j] + v[j] > cur_max){
        cur_max = values[i-j] + v[j];
      }
    }
    values[i] = cur_max;
  }
  return values[len];
}

// output the patterns of cutting that can generates the maximum value
// pred[] keeps track of the cut of each length
// pred[i] denotes the cut made when the rod length is i
int* rod_cutting_witness(int v[], int vsize, int len){
  int values[len+1];
  values[0] = 0;
  int cur_max;
  int *witness = new int [len+1];
  witness[0] = 0;
  for (int i=1; i<=len; i++){
    cur_max = -1;
    for (int j=1; j<=i; j++){
      if (j <vsize && values[i-j] + v[j] > cur_max){
        cur_max = values[i-j]+v[j];
        values[i] = cur_max;
        witness[i] = j;
      }
    }
  }
  int m = len;
  while (m != 0){
    if (witness[m] < m){
      cout<<"cut: "<<witness[m]<<" ";
      m = m - witness[m];
    }
    else {
      cout<<"cut: "<<witness[m]<<" ";
      break;
    }
  }
  cout<<endl;
  return witness;
}

int main(){
  int values1[] = {0,1,5,8,9,10,17,17,20};
  cout<<rod_cut(values1, sizeof(values1)/sizeof(values1[0]), 8)<<endl;
  cout<<rod_cut(values1, sizeof(values1)/sizeof(values1[0]), 20)<<endl;
  int *witness1_1 = rod_cutting_witness(values1, sizeof(values1)/sizeof(values1[0]), 8);
  int *witness1_2 = rod_cutting_witness(values1, sizeof(values1)/sizeof(values1[0]), 20);
  delete(witness1_1);
  delete(witness1_2);
}
