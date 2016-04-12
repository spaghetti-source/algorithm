// Minimum Coin Change problem (Dynamic Programming)
// Description:
//    Given an infinite supply of valued coins {c1, c2, c3, ... cn} and a value N
//    Output the minimum number of coins to make the change

// Algorithm:
//  Construct a table[][] and solve in a top bottom manner.
//  The witness pred[] array records the coins that are selected.
// define d[i] as the minimum number of coins to change i amout of money, cj as the coin value of jth coin
// then
//   d[i] = 0, when i=0
//   d[i] = min(d[i], d[i-cj]+1) for cj<=i, when i>=1

// Complexity:
// For m number of coins and n amount of money:
//  O(mn) for time
//  O(n) for space

#include <iostream>
#include <vector>
#include <cmath>
#include <climits>

using namespace std;

// find the minimum number of coin changes and print the witness
int coin_change(int coins[], int cn, int money){
  int table[money+1];
  table[0] = 0;
  int pred[money+1];
  for (int i=0; i<=money;i++){
    pred[i] = 0;
  }
  for (int j=1; j<=money;j++){
    table[j] = INT_MAX;
  }
  for (int i=1;i<=money;i++){
    int mini = table[i];
    for (int j=0; j<cn;j++){
      if (i >= coins[j]){
        mini = min(mini, table[i-coins[j]]+1);
        pred[i] = j;
      }
    }
    table[i] = mini;
  }
  int m = money;
  while (m != 0){
    cout<<"change coin: "<<coins[pred[m]]<<" ";
    m = m - coins[pred[m]];
  }
  cout<<endl;
  return table[money];
}

int main(){

  int coins1[] = {1, 5, 10, 25};
  cout<<coin_change(coins1, sizeof(coins1)/sizeof(coins1[0]), 20)<<" coins"<<endl;
  cout<<coin_change(coins1, sizeof(coins1)/sizeof(coins1[0]), 7)<<" coins"<<endl;

  int coins2[] = {1,5,7,15};
  cout<<coin_change(coins2, sizeof(coins2)/sizeof(coins2[0]), 8)<<" coins"<<endl;
  cout<<coin_change(coins2, sizeof(coins2)/sizeof(coins2[0]), 20)<<" coins"<<endl;
  return 0;
}
