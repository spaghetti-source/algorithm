/*
** Provide the value of n
** This program will print one of the solution of N-Queen problem if possible for given n.
** N-Queen problem is a great example of backtracking.
** where you chose a path and if the solution is not found you backtrack to choose a different path.
*/

#include<iostream>
using namespace std;
bool grid[10001][10001];
int n;
bool isAttacked(int x, int y)
{
    for(int i = 0; i < n; i++)
    {
        if(grid[i][y] || grid[x][i])
            return true;
    }
    for(int i = x-1, j = y-1; i >= 0; i--, j--)
    {
        if(grid[i][j])
            return true;
    }
    for(int i = x-1, j = y+1; i >= 0; i--, j++)
    {
        if(grid[i][j])
            return true;
    }
    return false;
}
bool nQueen(int row)
{
    if(row >= n)
        return true;
    for(int i = 0; i < n; i++)
    {
        if(!isAttacked(row, i))
        {
            grid[row][i] = true;
            if(nQueen(row + 1))
                return true;
            grid[row][i] = false;   //BACKTRACKING
        }
    }
    return false;
}
int main()
{
    cout<<"Enter number of queens : ";
    cin>>n;
    for(int i =0; i < n; i++)
        for(int j = 0; j < n; j++)
            grid[i][j] =false;
    if(nQueen(0))
    {
        cout<<"Solution :-\n";
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(grid[i][j])
                    cout<<"Q ";
                else
                    cout<<"- ";
            }
            cout<<endl;
        }
    }
    else
    {
        cout<<"No Answer.\n";
    }

    return 0;
}
