#include<iostream> 
using namespace std; 
  
// Returns value of Binomial Coefficient C(n, k) 
unsigned long int binomialCoeff(unsigned int n, unsigned int k) 
{ 
    unsigned long int res = 1; 
  
    // Since C(n, k) = C(n, n-k) 
    if (k > n - k) 
        k = n - k; 
  
    // Calculate value of [n*(n-1)*---*(n-k+1)] / [k*(k-1)*---*1] 
    for (int i = 0; i < k; ++i) 
    { 
        res *= (n - i); 
        res /= (i + 1); 
    } 
  
    return res; 
} 
  
// A Binomial coefficient based function to find nth catalan 
// number in O(n) time 
unsigned long int catalan(unsigned int n) 
{ 
    // Calculate value of 2nCn 
    unsigned long int c = binomialCoeff(2*n, n); 
  
    // return 2nCn/(n+1) 
    return c/(n+1); 
} 
  
// Driver program to test above functions 
int main() 
{ 
    for (int i = 0; i < 10; i++) 
        cout << catalan(i) << " "; 
    return 0; 
} 
