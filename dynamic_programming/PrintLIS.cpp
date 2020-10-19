void printLIS(vector<int>a, int n){
	vector<vector<int>>v(n);
	v[0].pb(a[0]);
	for(int i=1; i<n; i++){
		for(int j=0; j<i; j++){
			if(a[i]>a[j] && v[i].sz()<v[j].sz()+1)
				v[i]=v[j];
		}
		v[i].pb(a[i]);
	}
	int mx=0; int k;
	for(int i=0; i<n; i++)
		if(v[i].sz()>mx)
			mx=v[i].sz(), k=i;
	cout<<v[k].sz()<<nl;
	for(auto u:v[k])
		cout<<u<<" ";	
}
