
// A O(n) solution that uses  
// table fact[] to calculate  
// the Permutation Coefficient 
#include<bits/stdc++.h> 
  
// Returns value of Permutation 
// Coefficient P(n, k) 
int permutationCoeff(int n, int k) 
{ 
    int fact[n + 1]; 
  
    // base case 
    fact[0] = 1; 
  
    // Caculate value  
    // factorials up to n 
    for (int i = 1; i <= n; i++) 
        fact[i] = i * fact[i - 1]; 
  
    // P(n,k) = n! / (n - k)! 
    return fact[n] / fact[n - k]; 
} 
  
// Driver Code 
int main() 
{ 
    int n = 10, k = 2; 
    printf ("Value of P(%d, %d) is %d ", 
             n, k, permutationCoeff(n, k) ); 
    return 0; 
} 
