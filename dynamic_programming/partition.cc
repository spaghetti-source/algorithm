// Partition Problem (Dynamic Programming)

// Description:
//  Given an integer array a, determine whether a can be partitioned into two subsets
//  s1 and s2, such that sum(s1) = sum(s2)

// Example:
//  input: arr[] = {4,1,2,3}
//  output: True
//  The array can be partitioned to {4,1} and {2, 3}

#include <iostream>

using namespace std;

bool partition(int arr[], int n){
  int sum = 0;
  for (int i=0; i<n ;i++){
    sum += arr[i];
  }
  if (sum % 2 != 0){
    return false;
  }
  sum /= 2;
  int matrix[n+1][sum+1];
  for (int i=0; i<=n; i++){
    for (int j=0; j<=sum; j++){
      if (j==0){
        matrix[i][j] = true;
      }
      else if (i==0){
        matrix[i][j] = false;
      }
      else {
        matrix[i][j] = matrix[i-1][j];
        if (j >= arr[i-1]){
          matrix[i][j] = matrix[i][j] || matrix[i-1][j-arr[i-1]];
        }
      }
    }
  }
  return matrix[n][sum];
}
int main(){
  int arr1[] = {3,1,1,2,2,1};
  int arr2[] = {4,2,1,3,4};
  int arr3[] = {2,4,5,5};
  cout<<partition(arr1, sizeof(arr1)/sizeof(arr1[0]))<<endl;
  cout<<partition(arr2, sizeof(arr2)/sizeof(arr2[0]))<<endl;
  cout<<partition(arr3, sizeof(arr3)/sizeof(arr3[0]))<<endl;
  return 0;
}
