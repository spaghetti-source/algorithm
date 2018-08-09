#include <bits/stdc++.h>

using namespace std;

double trapezoid(double x1, double x2, double y1, double y2) {
	return (y2+y1)/2 * abs(x2-x1);
}

double auc(vector<int> test, vector<double> pred) {
	int n = test.size();
	assert(n == pred.size());

	vector<int> idx(n);
	for (int i = 0; i < n; ++i) idx[i] = i;
	sort(idx.begin(), idx.end(), [&](int i, int j) { return pred[i] > pred[j]; });

	double a = 0.0;
	double fp = 0, tp = 0, fp_prev = 0, tp_prev = 0;
	double prev_score = -1.0/0.0;
	for (int i: idx) {
		if (pred[i] != prev_score) {
			a += trapezoid(fp, fp_prev, tp, tp_prev);
			prev_score = pred[i];
			fp_prev = fp;
			tp_prev = tp;
		}
		if (test[i] == 1) {
			tp += 1;
		} else {
			fp += 1;
		}
	}
	a += trapezoid(fp, fp_prev, tp, tp_prev);
	return a / (tp * fp);
}
int main() {
	vector<int> test = {0, 1, 0, 1, 1};
	vector<double> pred = {0.2, 0.3, 0.4, 0.5, 0.6};
	cout << auc(test, pred) << endl;
}
