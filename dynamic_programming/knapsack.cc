#include<iostream>
using namespace std;

int main()
{
    int n,W,i,j,count=0;
    cout<<"Enter the weight capacity of the bag : ";
    cin>>W;
    cout<<"Enter the number of items :";
    cin>>n;
    int v[n],w[n],T[n+1][W+1];
    cout<<"Enter weight and value of each item respectively :\n";
    for(i=0;i<n;i++)
    {
        cout<<i+1<<" ";
        cin>>w[i]>>v[i];
    }

    for(j=0;j<=W;j++)
        T[0][j]=0;
    for(i=0;i<=n;i++)
        T[i][0]=0;
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=W;j++)
        {
            if(w[i-1]<=j)
                if(T[i-1][j]<v[i-1]+T[i-1][j-w[i-1]])
                    T[i][j]=v[i-1]+T[i-1][j-w[i-1]];
                else
                    T[i][j]=T[i-1][j];
            else
                T[i][j]=T[i-1][j];
        }
    }

    bool req[n+1];

    for(i=0;i<n+1;i++)
        req[i]=0;

    for(j=W-1;j>0;j--)
        for(i=n-1;i>0;i--)
            if(T[i+1][j]!=T[i][j])
                req[i]=1;
    cout<<"The maximum value is : "<<T[n][W]<<"\nThe items to be put into knapsack : ";
    for(i=1;i<=n;i++)
        if(req[i])
            cout<<i<<"\t";

    return 0;

}
