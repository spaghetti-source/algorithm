#include<bits/stdc++.h>
#define fr(siz,i) for(int i=0;i<siz;i++)
#define frr(siz,i,a) for(int i=a;i<siz;i++)
#define ff(siz) for(int i=0;i<siz;i++)
#define ll long long
#define pb push_back
#define pii pair<int,int>
#define pll pair< ll , ll >
#define vi vector<int>
#define vvi vector< vi >
#define vl vector<ll>
#define vvl vector< vl >
const int maxn=(int)(2e5+5);
const ll mod=(ll)(1e9+7);
//ios_base::sync_with_stdio(0);cin.tie(0);
using namespace std;


struct trie_node{
	char data;
	bool isEnd;
	map<char,trie_node*> children;
};

trie_node * newNode(char data,bool isEnd) {
	trie_node* temp = new trie_node;
	temp->data = data;
	temp->isEnd = isEnd;
	return temp;
}

bool insert(trie_node* head, string data) {
	trie_node* currentNode = head;
	bool wasThere = true;
	for(int i=0;i<data.length();i++) {
		if(currentNode->children[data[i] - 'a'] == NULL) {
			bool isEnd = ( i == data.length()-1 );
			trie_node* temp = newNode(data[i],isEnd);
			currentNode->children[data[i] - 'a'] = temp;
			wasThere = false;
		}
		currentNode = currentNode->children[data[i] - 'a'];
	}
	return wasThere;
}

bool find(trie_node* head, string data, bool prefix) {
	trie_node* currentNode = head;
	bool found = true;
	for(int i=0;i<data.length();i++) {
		if(currentNode->children[data[i] - 'a'] == NULL) {
			return false;
		}
		currentNode = currentNode->children[data[i] - 'a'];
	}
	return prefix ? true : currentNode->isEnd ;
}

int main() {

	trie_node * head = newNode('*',true);
    cout<<insert(head,"sasta_achar");
	cout<<insert(head,"sasta_achar");
	cout<<insert(head,"sasta_pickle");
	cout<<"string"<<find(head,"don",false);
 	cout<<"\nprefix don"<<find(head,"don",true);
	return 0;
}
