#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include "ParabolMethod.h"

#define DELTA 0.05

//============================================================PROTECTED=========================================================

void ParabolMethod::calculateR(const shared_ptr<Segment>& s) {
	switch (mode) {
	case 1:
		calculateRAnalytics(s);
		break;
	case 2:
		calculateRSections(s);
		break;
	default:
		break;
	}
}


void ParabolMethod::calculateRAnalytics(const shared_ptr<Segment>& s) {

	if (accVal == true) {
		auto p1 = getNegAreaOfConstraints(s);

		pair<double, double> f_ = pair<double, double>(INFINITY, INFINITY);
		// �������� �� ����������� ����������� ���������
		for (int i = 0; i < p1.size(); i++) {
			pair<double, double> newf_ = findMin(s, p1[i]);
			if (newf_.second <= f_.second) { // ������ �������������
				f_.first = newf_.first;
				f_.second = newf_.second;
			}
		}
		s->R = f_.second;
		s->argmin = f_.first;
	}
	else {
		auto p1 = getNegAreaOfConstraints(s, currDiscrep);
		pair<double, double> f_ = pair<double, double>(INFINITY, INFINITY);
		// �������� �� ����������� ����������� ���������
		for (int i = 0; i < p1.size(); i++) {
			pair<double, double> newf_ = fibonacciMethod(s, p1[i]);
			if (newf_.second <= f_.second) { // ������ �������������
				f_.first = newf_.first;
				f_.second = newf_.second;
			}
		}
		s->R = f_.second;
		s->argmin = f_.first;
	}
}

 //���������� �������������� ��������� �� ����� �����
void ParabolMethod::calculateRSections(const shared_ptr<Segment>& s) {
	const int pointCoint = 300; // ����� ����� "��������" �� ���������

	vector<double> minInPoint; // ������� ������� - �������� �� ����������� �������� ������ ������ � 1 �����
	minInPoint.resize(pointCoint);

	double alpha = s->length / (pointCoint + 1);
	double point;

	if (accVal == true) {

		for (int i = 0; i < pointCoint; i++) {
			point = s->p1->x + alpha * (i + 1);

			bool Need = true;
			// ���� ���� �� ���� ����������� ������ ���� �� ��������� ����� ������
			for (int p = 1; p < constraintsCount + 1; ++p) {
				if (point < s->xt1[p]) {
					//���� ��� �������� ������� ����, �����
					minInPoint[i] = s->ca[p] - (s->L[p] / 2) * (point - s->wa[p]) * (point - s->wa[p]);
				}
				else if (point > s->xt2[p]) {
					//���� ��� �������� ������� ����, ������
					minInPoint[i] = s->cb[p] - (s->L[p] / 2) * (point - s->wb[p]) * (point - s->wb[p]);
				}
				else if (point <= s->xt2[p] && point >= s->xt1[p]) {
					//���� ��� �������� ������� �����, � ������
					minInPoint[i] = s->c[p] + (s->L[p] / 2) * (point - s->w[p]) * (point - s->w[p]);
				}
				if (minInPoint[i] >= 0) {
					Need = false;
					minInPoint[i] = INFINITY;
					break;
				}
			}
			if (Need != false) {
				if (point < s->xt1[0]) {
					//���� ��� �������� ������� ����, �����
					minInPoint[i] = s->ca[0] - (s->L[0] / 2) * (point - s->wa[0]) * (point - s->wa[0]);
				}
				else if (point > s->xt2[0]) {
					//���� ��� �������� ������� ����, ������
					minInPoint[i] = s->cb[0] - (s->L[0] / 2) * (point - s->wb[0]) * (point - s->wb[0]);
				}
				else if (point <= s->xt2[0] && point >= s->xt1[0]) {
					//���� ��� �������� ������� �����, � ������
					minInPoint[i] = s->c[0] + (s->L[0] / 2) * (point - s->w[0]) * (point - s->w[0]);
				}
			}
		}

		int argminpoint = 0;
		//������ ���������� maxInPoint �� ���� ������ ��������
		//������� ����� ���
		double min = minInPoint[0];
		for (int i = 1; i < pointCoint; i++) {
			if (min > minInPoint[i]) {
				min = minInPoint[i];
				argminpoint = i;
			}
		}
		//�������������� ��������� - ��������� �������
		s->R = min;
		s->argmin = s->p1->x + alpha * (argminpoint + 1);
	}
	else {

		for (int i = 0; i < pointCoint; i++) {
			point = s->p1->x + alpha * (i + 1);

			bool Need = true;
			double MaxConstraint = -INFINITY;
			// ���� ���� �� ���� ����������� ������ ���� �� ��������� ����� ������
			for (int p = 1; p < constraintsCount + 1; ++p) {
				if (point < s->xt1[p]) {
					//���� ��� �������� ������� ����, �����
					minInPoint[i] = s->ca[p] - (s->L[p] / 2) * (point - s->wa[p]) * (point - s->wa[p]);
				}
				else if (point > s->xt2[p]) {
					//���� ��� �������� ������� ����, ������
					minInPoint[i] = s->cb[p] - (s->L[p] / 2) * (point - s->wb[p]) * (point - s->wb[p]);
				}
				else if (point <= s->xt2[p] && point >= s->xt1[p]) {
					//���� ��� �������� ������� �����, � ������
					minInPoint[i] = s->c[p] + (s->L[p] / 2) * (point - s->w[p]) * (point - s->w[p]);
				}
				if (minInPoint[i] >= currDiscrep) {
					Need = false;
					minInPoint[i] = INFINITY;
					break;
				}
				// �������� ������� ���������
				MaxConstraint = max(minInPoint[i], MaxConstraint);
			}
			if (Need != false) {
				minInPoint[i] = MaxConstraint;
			}
		}

		int argminpoint = 0;
		//������ �������� ������� ��������� �� ���� ������ ��������
		//������� ����� ���
		double min = minInPoint[0];
		for (int i = 1; i < pointCoint; i++) {
			if (min > minInPoint[i]) {
				min = minInPoint[i];
				argminpoint = i;
			}
		}
		//�������������� ��������� - ��������� �������
		s->R = min;
		currDiscrep = min;
		s->argmin = s->p1->x + alpha * (argminpoint + 1);
	}
}

void ParabolMethod::calculateLowerLimit(const shared_ptr<Segment>& s) {
	for (size_t i = 0; i < constraintsCount + 1; i++) {
		double x = s->p1->x;
		double fx = s->p1->f[i]; //�������� i-�� ������� � ����� �
		double fxD = s->p1->fd[i]; //�������� ����������� i-�� ������� � ����� �
		s->wa[i] = (fxD + x * s->L[i]) / s->L[i];
		s->ca[i] = s->L[i] / 2 * (-1) * pow(s->wa[i] - s->p1->x, 2) + fxD * (s->wa[i] - s->p1->x) + fx;

		x = s->p2->x;
		fx = s->p2->f[i];
		fxD = s->p2->fd[i];
		s->wb[i] = (fxD + x * s->L[i]) / s->L[i];
		s->cb[i] = s->L[i] / 2 * (-1) * pow(s->wb[i] - s->p2->x, 2) + fxD * (s->wb[i] - s->p2->x) + fx;

		s->xt1[i] = s->wa[i] - (s->wa[i] - s->wb[i]) / 4 + (s->cb[i] - s->ca[i]) / (s->L[i] * (s->wa[i] - s->wb[i]));
		s->xt2[i] = s->wb[i] + (s->wa[i] - s->wb[i]) / 4 + (s->cb[i] - s->ca[i]) / (s->L[i] * (s->wa[i] - s->wb[i]));
		//s->w[i] = (s->ca[i] - s->cb[i]) / (2 * s->L[i] * (s->xt2[i] - s->xt1[i])) + (s->xt2[i] + s->xt1[i]) / 2;
		s->w[i] = (s->wb[i] + s->wa[i]) / 2 - 2 * ((s->cb[i] - s->ca[i]) / ((s->wb[i] - s->wa[i]) * s->L[i]));
		//s->c[i] = (s->ca[i] + s->cb[i] - s->L[i] * pow((s->xt1[i] - s->wa[i]), 2) - s->L[i] * pow((s->xt2[i] - s->wb[i]), 2) / 2);
		s->c[i] = (s->ca[i] + s->cb[i]) / 2 -
			s->L[i] / 16 * pow(s->wb[i] - s->wa[i], 2) -
			pow(s->cb[i] - s->ca[i], 2) / (pow(s->wb[i] - s->wa[i], 2) * s->L[i]);
	}
}

double ParabolMethod::findLLoc(const shared_ptr<Segment>& s, const int i) const {
	//calculateM
	double fd1 = s->p1->fd[i];
	double fd2 = s->p2->fd[i];
	double f1 = s->p1->f[i];
	double f2 = s->p2->f[i];
	double value, value1;
	value1 = ((f2 - f1) - (s->p2->x - s->p1->x) * ((fd2 + fd1) / 2)) / pow((s->p2->x - s->p1->x), 2);
	value = 2 * abs(value1) + sqrt(4 * pow(value1, 2) + pow(((fd2 - fd1) / (s->p2->x - s->p1->x)), 2));
	return value;
}

//============================================================PRIVATE===========================================================

vector<double> ParabolMethod::findRoots(const double a, const double b, const double c) const {
	vector<double> x;
	x.reserve(2);

	double d = pow(b, 2) - 4 * a * c;
	if (d > 0) {
		double x1 = (-b + sqrt(d)) / (2 * a);
		double x2 = (-b - sqrt(d)) / (2 * a);
		x.push_back(x1 < x2 ? x1 : x2);
		x.push_back(x1 < x2 ? x2 : x1);
	}
	else if (d == 0) {
		x.push_back(-b / (2 * a));
	}
	return x;
}

vector<pair<double, double>> ParabolMethod::getNegAreaOfConstraints(const shared_ptr<Segment>& s) const {
	vector<pair<double, double>> p1{make_pair(s->p1->x, s->p2->x)};
	vector<double> rootsA, rootsB, roots;
	vector<pair<double, double>> v1, v2, v3;
	// ��� ������ ������� ����������� ���� ����� �� ��������� � ���������� � ���������� ������� �������
	for (int i = 1; i < constraintsCount + 1; i++) {
		// ��� ������ ������ ��������� �� ���� ������� �� ���� ������� ������� ����� ������ ��������
		rootsA = findRoots(-s->L[i] / 2, s->L[i] * s->wa[i], s->ca[i] - (s->L[i] / 2) * pow(s->wa[i], 2));
		rootsB = findRoots(-s->L[i] / 2, s->L[i] * s->wb[i], s->cb[i] - (s->L[i] / 2) * pow(s->wb[i], 2));
		roots = findRoots(s->L[i] / 2, -s->L[i] * s->w[i], s->c[i] + (s->L[i] / 2) * pow(s->w[i], 2)); //��������� ������� C + M/2(x - V)^2
		// ������� �� ����� ��������� ������ �� ������� ��������� ������������� ��������
		v1 = validSetOfNegPar(s->p1->x, s->xt1[i], rootsA);
		v2 = validSetOfNegPar(s->xt2[i], s->p2->x, rootsB);
		v3 = validSetOfPosPar(s->xt1[i], s->xt2[i], roots);
		// ���������� ��������� ��� �������� ��������� ������������� �������� � ���� ���������
		vector<pair<double, double>> p2 = combine(v1, v2, v3);
		p1 = intersectUnions(p1, p2);
	} // �������� ��������� ��� ��� ������ ������ ������� ��������� ������������� �������� - ������ ����������� ���������, ��� ����� ������ �������
	return p1;
}

vector<pair<double, double>> ParabolMethod::getNegAreaOfConstraints(const shared_ptr<Segment>& s, double param) const {
	vector<pair<double, double>> p1{ make_pair(s->p1->x, s->p2->x) };
	vector<double> rootsA, rootsB, roots;
	vector<pair<double, double>> v1, v2, v3;
	for (int i = 1; i < constraintsCount + 1; i++) {
		// ��� ������ ������ ��������� �� ���� ������� �� ���� ������� ������� ����� ������ �������� � ������ ������ �������
		rootsA = findRoots(-s->L[i] / 2, s->L[i] * s->wa[i], s->ca[i] - (s->L[i] / 2) * pow(s->wa[i], 2) - param);
		rootsB = findRoots(-s->L[i] / 2, s->L[i] * s->wb[i], s->cb[i] - (s->L[i] / 2) * pow(s->wb[i], 2) - param);
		roots = findRoots(s->L[i] / 2, -s->L[i] * s->w[i], s->c[i] + (s->L[i] / 2) * pow(s->w[i], 2) - param); //��������� ������� C + M/2(x - V)^2
		// ������� �� ����� ��������� ������ �� ������� ��������� �������� ������ ������ �������
		v1 = validSetOfNegPar(s->p1->x, s->xt1[i], rootsA);
		v2 = validSetOfNegPar(s->xt2[i], s->p2->x, rootsB);
		v3 = validSetOfPosPar(s->xt1[i], s->xt2[i], roots);
		// ���������� ��������� ��� �������� ��������� �������� ������ ������ ������� � ���� ���������
		vector<pair<double, double>> p2 = combine(v1, v2, v3);
		p1 = intersectUnions(p1, p2);
	} // �������� ��������� ��� ��� ������ ������ ������� ��������� �������� ������ �������� ������ �������
	return p1;
}



vector<pair<double, double>> ParabolMethod::validSetOfNegPar(const double a, const double b, const vector<double> roots) const {
	vector<pair<double, double>> res; // ��� ���������� ���, ��� ������ roots[0] � ���, ��� ������ roots[1], ������ ���� 6 �������
	res.reserve(2);
	if (b < a) {
		if (roots.empty()) {
			res.push_back(pair<double, double>(b, a));
			return res; 
		}
		else {
			return res;
		}
	}
	if (roots.empty() || roots.size() == 1) {
		res.push_back(pair<double, double>(a, b));
		return res;
	}
	if (roots[0] < a && roots[1] > b) {
		return res;
	}
	if (roots[0] > a && roots[0] < b && roots[1] < b) { // ������, ����� ��� ����� �� p1 �� x'
		res.push_back(pair<double, double>(a, roots[0]));
		res.push_back(pair<double, double>(roots[1], b));
	}
	if (roots[0] < a && roots[1] > a && roots[1] < b) { // roots[0] ������� �� �������, roots[1] � ���������
		res.push_back(pair<double, double>(roots[1], b));
	}
	if (roots[0] > a && roots[1] > b && roots[0] < b) { // roots[1] ������� �� �������, roots[0] � ���������
		res.push_back(pair<double, double>(a, roots[0]));
	}
	if (roots[0] < a && roots[1] <= a) {
		res.push_back(pair<double, double>(a, b)); // roots[0] ������� �� �������, roots[1] �� ����� �������
	}
	if (roots[0] >= b && roots[1] > b) { // roots[1] ������� �� �������, roots[0] �� �������
		res.push_back(pair<double, double>(a, b));
	}
	return res;
}

vector<pair<double, double>> ParabolMethod::validSetOfPosPar(const double a, const double b, const vector<double> roots) const {
	vector<pair<double, double>> res; //  ������� x' � x'', ���������� ���, ��� ������ 0, ������� ������� roots[0] � ������� roots[1]
	double first = a;
	double second = b;
	res.reserve(2);
	if (b < a) {
		if (roots.empty()) {
			return res;
		}
		else {
			first = b;
			second = a; 
		}
	}
	if (roots.empty() ) {
		return res;
	}
	if (roots.size() == 1) {
		if (roots[0] >= first && roots[0] <= second) {
			res.push_back(pair<double, double>(roots[0], roots[0]));
		}
		return res;
	}
	int count = 0;
	if (roots[0] <= first && roots[1] >= second) { // roots[i] �� �������� ��� ������, ����� ���� ��������
		res.push_back(pair<double, double>(first, second));
	}
	if (roots[0] < first && roots[1] < second && roots[1] > first) { // roots[0] �� ��������, roots[1] � ���������
		res.push_back(pair<double, double>(first, roots[1]));
	}
	if (first < roots[0] && roots[0] < second && roots[1] > second) { // roots[1] �� ��������, roots[0] � ���������
		res.push_back(pair<double, double>(roots[0], second));
	}
	if (roots[0] >= first && roots[1] <= second) { // roots[0] roots[1] ������ ���������
		res.push_back(pair<double, double>(roots[0], roots[1]));
	}
	return res; // ����� ������� ������ ����������� ���� �������� ����� !!
}

vector<pair<double, double>> ParabolMethod::combine(const vector<pair<double, double>> v1, const vector<pair<double, double>> v2, const vector<pair<double, double>> v3) const {
	vector<pair<double, double>> tmp;
	vector<pair<double, double>> result;
	pair<double, double> s;

	for (int i = 0; i < v1.size(); i++)
	{
		tmp.push_back(v1[i]);
	}
	for (int i = 0; i < v2.size(); i++)
	{
		tmp.push_back(v2[i]);
	}
	for (int i = 0; i < v3.size(); i++)
	{
		tmp.push_back(v3[i]);
	}
	if (tmp.empty()) {
		return tmp;
	}

	bool flag = true;
	int lastElem = 0;
	while (!tmp.empty())
	{
		result.push_back(tmp[0]);
		tmp.erase(tmp.begin());
		lastElem = result.size() - 1;
		flag = true;
		while (flag) {
			flag = false;
			for (int j = 0; j < tmp.size(); j++)
			{
				if (result[lastElem].second == tmp[j].first)
				{
					result[lastElem] = pair<double, double>(result[lastElem].first, tmp[j].second);
					tmp.erase(tmp.begin() + j);
					flag = true;
					break;
				}
				else if (tmp[j].second == result[lastElem].first)
				{
					result[lastElem] = pair<double, double>(tmp[j].first, result[lastElem].second);
					tmp.erase(tmp.begin() + j);
					flag = true;
					break;
				}
			}
		}
	}


	return result;
}

pair<double, double> ParabolMethod::minPos(const double x, const double y, const double a,
	                                       const double b, const double c, const double arg_min) const
{
	double argMin, min = INFINITY;
	double fx = a * pow(x, 2) + b * x + c;
	if (fx < min) {
		min = fx;
		argMin = x;
	}
	fx = a * pow(y, 2) + b * y + c;
	if (fx < min) {
		min = fx;
		argMin = y;
	}
	if (arg_min >= x && arg_min <= y) {

		fx = a * pow(arg_min, 2) + b * arg_min + c;
		if (fx < min) {
			min = fx;
			argMin = arg_min;
		}
	}
	return pair<double, double>(argMin, min);
}

pair<double, double> ParabolMethod::minNeg(const double x, const double y, const double a, const double b, const double c) const
{
	double argMin, min = INFINITY;
	double fx = a * pow(x, 2) + b * x + c;
	if (fx < min) {
		min = fx;
		argMin = x;
	}
	fx = a * pow(y, 2) + b * y + c;
	if (fx < min) {
		min = fx;
		argMin = y;
	}
	return pair<double, double>(argMin, min);
}


vector<pair<double, double>> ParabolMethod::intersectUnions(const vector<pair<double, double>>x, const vector<pair<double, double>>y) const {
	vector<pair<double, double>> res;
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < y.size(); j++) {
			pair<double, double> s = intersectSegments(x[i], y[j]);
			if (s.first != -INFINITY && s.second != INFINITY) {
				res.push_back(s);
			}
		}
	}
	return res;
}

pair<double, double> ParabolMethod::findMin(const shared_ptr<Segment> s, const pair<double, double> p) const {
	pair<double, double> min(INFINITY, INFINITY);
	pair<double, double> min1;
	pair<double, double> min2;
	pair<double, double> min3;
	if (p.second <= s->xt1[0] && p.first <= s->xt1[0]) { // 1 ��������
		return minNeg(p.first, p.second, -s->L[0] / 2, s->L[0] * s->wa[0], s->ca[0] - (s->L[0] / 2) * pow(s->wa[0], 2));
	}
	if (p.first >= s->xt1[0] && p.second <= s->xt2[0]) { // 2 �������� 
		return minPos(p.first, p.second, s->L[0] / 2, -s->L[0] * s->w[0], s->c[0] + (s->L[0] / 2) * pow(s->w[0], 2), s->w[0]);
	}
	if (p.first >= s->xt2[0] && p.second >= s->xt2[0]) { // 3 ��������
		return minNeg(p.first, p.second, -s->L[0] / 2, s->L[0] * s->wb[0], s->cb[0] - (s->L[0] / 2) * pow(s->wb[0], 2));
	}
	if (p.first < s->xt1[0] && p.second > s->xt1[0] && p.second < s->xt2[0]) { // ����� ����� 1 � 2 ��������
		min1 = minNeg(p.first, s->xt1[0], -s->L[0] / 2, s->L[0] * s->wa[0], s->ca[0] - (s->L[0] / 2) * pow(s->wa[0], 2));
		min2 = minPos(s->xt1[0], p.second, s->L[0] / 2, -s->L[0] * s->w[0], s->c[0] + (s->L[0] / 2) * pow(s->w[0], 2), s->w[0]);
		return min1.second < min2.second ? min1 : min2;
	}
	if (p.first > s->xt1[0] && p.first < s->xt2[0] && p.second > s->xt2[0]) { // ����� 2 � 3
		min3 = minNeg(s->xt2[0], p.second, -s->L[0] / 2, s->L[0] * s->wb[0], s->cb[0] - (s->L[0] / 2) * pow(s->wb[0], 2));
		min2 = minPos(p.first, s->xt2[0], s->L[0] / 2, -s->L[0] * s->w[0], s->c[0] + (s->L[0] / 2) * pow(s->w[0], 2), s->w[0]);
		return min2.second < min3.second ? min2 : min3;
	}
	if (p.first < s->xt1[0] && p.second > s->xt2[0]) { // ����� ����� ���� ����
		min1 = minNeg(p.first, s->xt1[0], -s->L[0] / 2, s->L[0] * s->wa[0], s->ca[0] - (s->L[0] / 2) * pow(s->wa[0], 2));
		min2 = minPos(s->xt1[0], s->xt2[0], s->L[0] / 2, -s->L[0] * s->w[0], s->c[0] + (s->L[0] / 2) * pow(s->w[0], 2), s->w[0]);
		min3 = minNeg(s->xt2[0], p.second, -s->L[0] / 2, s->L[0] * s->wb[0], s->cb[0] - (s->L[0] / 2) * pow(s->wb[0], 2));
		min = min1.second < min2.second ? min1 : min2;
		min = min.second < min3.second ? min : min3;
		return min;
	}
	return min;
}

// TODO remove
void ParabolMethod::printViewData()
{
	int lineZero = 0;
	std::ofstream out("viewDataParabol.txt", ios_base::trunc);
	int pointsCount = 200;
	vector<shared_ptr<Segment>>::iterator IT = segments.begin();
	double h = 1.0 / (pointsCount - 1u), x;
	for (int i = 0; i < pointsCount; ++i)
	{
		x = i * h;
		if (x > (*IT)->p2->x)
		{
			++IT;
		}
		//���������� ������� ��������� �������:
		TProblem& P = (*t);
		double maxFunc = P(x, 0);
		double maxConstr = P(x, 1);

	//���������� ������� ��������� �������� �����������  
		double maxFuncMinor = calculateMinorant((*IT), x);
	// ���������� ������ ������ ������� ������� ��� ��������� ��� ������ ����� ����������� <=0
		double InPoint;
		bool Need = true;
				// ���� ���� �� ���� ����������� ������ ���� �� ��������� ����� ������
				for (int p = 1; p < constraintsCount + 1; ++p) {
					if (x < (*IT)->xt1[p]) {
						//���� ��� �������� ������� ����, �����
						InPoint = (*IT)->ca[p] - ((*IT)->L[p] / 2) * (x - (*IT)->wa[p]) * (x - (*IT)->wa[p]);
					}
					else if (x > (*IT)->xt2[p]) {
						//���� ��� �������� ������� ����, ������
						InPoint = (*IT)->cb[p] - ((*IT)->L[p] / 2) * (x - (*IT)->wb[p]) * (x - (*IT)->wb[p]);
					}
					else if (x <= (*IT)->xt2[p] && x >= (*IT)->xt1[p]) {
						//���� ��� �������� ������� �����, � ������
						InPoint = (*IT)->c[p] + ((*IT)->L[p] / 2) * (x - (*IT)->w[p]) * (x - (*IT)->w[p]);
					}
					if (InPoint >= 0) {
						Need = false;
						InPoint = 0;
						break;
					}
				}
				if (Need != false) {
					if (x < (*IT)->xt1[0]) {
						//���� ��� �������� ������� ����, �����
						InPoint = (*IT)->ca[0] - ((*IT)->L[0] / 2) * (x - (*IT)->wa[0]) * (x - (*IT)->wa[0]);
					}
					else if (x > (*IT)->xt2[0]) {
						//���� ��� �������� ������� ����, ������
						InPoint = (*IT)->cb[0] - ((*IT)->L[0] / 2) * (x - (*IT)->wb[0]) * (x - (*IT)->wb[0]);
					}
					else if (x <= (*IT)->xt2[0] && x >= (*IT)->xt1[0]) {
						//���� ��� �������� ������� �����, � ������
						InPoint = (*IT)->c[0] + ((*IT)->L[0] / 2) * (x - (*IT)->w[0]) * (x - (*IT)->w[0]);
					}
				}

		out << x << " " << lineZero << " " << maxFunc << " " << maxConstr << " " << maxFuncMinor << " " << InPoint << std::endl;

	}
	out.close();
}

// TODO remove
double ParabolMethod::calculateMinorant(const shared_ptr<Segment>& s, double x)
{
	double result = -INFINITY; // ���������� ������ ������ ����� �� ������� �����������
	double maximum = -INFINITY; //�������� ������� ��������� ������ � ����� �����
		//�� ������ ������� �����������
	for (int p = 1; p < constraintsCount + 1; ++p) {
		double tmp1 = s->ca[p] - (s->L[p] / 2) * (s->xt1[p] - s->wa[p]) * (s->xt1[p] - s->wa[p]);
		double tmp2 = s->c[p] + (s->L[p] / 2) * (s->xt1[p] - s->w[p]) * (s->xt1[p] - s->w[p]);

		if (x < s->xt1[p]) {
			//���� ��� �������� ������� ����, �����
			result = s->ca[p] - (s->L[p] / 2) * (x - s->wa[p]) * (x - s->wa[p]);
		}
		else if (x > s->xt2[p]) {
			//���� ��� �������� ������� ����, ������
			result = s->cb[p] - (s->L[p] / 2) * (x - s->wb[p]) * (x - s->wb[p]);
		}
		else if (x <= s->xt2[p] && x >= s->xt1[p]) {
			//���� ��� �������� ������� �����, � ������
			result = s->c[p] + (s->L[p] / 2) * (x - s->w[p]) * (x - s->w[p]);
		}
		if (p != 0) {
			maximum = max(maximum, result);
		}
		else {
			maximum = max(maximum, result);
		}
	}
	return maximum;
}

void ParabolMethod::generateFibonacci(int n) {
	Fibonacci.push_back(1);
	Fibonacci.push_back(1);
	for (int i = 2; i < n; i++) {
		Fibonacci.push_back(Fibonacci[i - 1] + Fibonacci[i - 2]);
	}
}

pair<double, double> ParabolMethod::getEnvelope(const shared_ptr<Segment>& s, double x1, double x2) {
	double xlMax = -INFINITY, xrMax = -INFINITY;
	for (int i = 1; i < constraintsCount + 1; i++) {
		if (x2 <= s->xt1[i] && x1 <= s->xt1[i]) { // 1 ��������
			xlMax = s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) > xlMax ?
				s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) : xlMax;

			xrMax = s->ca[i] - (s->L[i] / 2) * (x2 - s->wa[i]) * (x2 - s->wa[i]) > xrMax ?
				s->ca[i] - (s->L[i] / 2) * (x2 - s->wa[i]) * (x2 - s->wa[i]) : xrMax;
		}
		if (x1 >= s->xt1[i] && x2 <= s->xt2[i]) { // 2 �������� 
			xlMax = s->c[i] - (s->L[i] / 2) * (x1 - s->w[i]) * (x1 - s->w[i]) > xlMax ?
				s->c[i] - (s->L[i] / 2) * (x1 - s->w[i]) * (x1 - s->w[i]) : xlMax;

			xrMax = s->c[i] - (s->L[i] / 2) * (x2 - s->w[i]) * (x2 - s->w[i]) > xrMax ?
				s->c[i] - (s->L[i] / 2) * (x2 - s->w[i]) * (x2 - s->w[i]) : xrMax;
		}
		if (x1 >= s->xt2[i] && x2 >= s->xt2[i]) { // 3 ��������
			xlMax = s->cb[i] - (s->L[i] / 2) * (x1 - s->wb[i]) * (x1 - s->wb[i]) > xlMax ?
				s->cb[i] - (s->L[i] / 2) * (x1 - s->wb[i]) * (x1 - s->wb[i]) : xlMax;

			xrMax = s->c[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) > xrMax ?
				s->cb[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) : xrMax;
		}
		if (x1 < s->xt1[i] && x2 > s->xt1[i] && x2 < s->xt2[i]) { // ����� ����� 1 � 2 ��������
			xlMax = s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) > xlMax ?
				s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) : xlMax;

			xrMax = s->c[i] - (s->L[i] / 2) * (x2 - s->w[i]) * (x2 - s->w[i]) > xrMax ?
				s->c[i] - (s->L[i] / 2) * (x2 - s->w[i]) * (x2 - s->w[i]) : xrMax;
		}
		if (x1 > s->xt1[i] && x1 < s->xt2[i] && x2 > s->xt2[i]) { // ����� 2 � 3
			xlMax = s->c[i] - (s->L[i] / 2) * (x1 - s->w[i]) * (x1 - s->w[i]) > xlMax ?
				s->c[i] - (s->L[i] / 2) * (x1 - s->w[i]) * (x1 - s->w[i]) : xlMax;

			xrMax = s->c[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) > xrMax ?
				s->cb[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) : xrMax;
		}
		if (x1 < s->xt1[i] && x2 > s->xt2[i]) { // ����� ����� ���� ����
			xlMax = s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) > xlMax ?
				s->ca[i] - (s->L[i] / 2) * (x1 - s->wa[i]) * (x1 - s->wa[i]) : xlMax;

			xrMax = s->c[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) > xrMax ?
				s->cb[i] - (s->L[i] / 2) * (x2 - s->wb[i]) * (x2 - s->wb[i]) : xrMax;
		}
	}
	return pair<double, double> {xlMax, xrMax};
}


pair<double, double> ParabolMethod::fibonacciMethod(const shared_ptr<Segment>& s, pair<double, double> seg) {
	double left = seg.first;
	double right = seg.second;
	double x1, x2;
	int n = ceil(log2((right - left) / eps)); // ���������� �����
	double deltaK = right - left;

	double min = INFINITY, argmin = INFINITY;
	pair<double, double> splice; // �������� ������������ ������ ������ �-� �����������

	// ����� ���������
	for (int j = 0; j < n; j++) 
	{
		deltaK = right - left;
		x1 = left + deltaK * (Fibonacci[n - 1 - j] / Fibonacci[n - j]);
		x2 = right - deltaK * (Fibonacci[n - 1 - j] / Fibonacci[n - j]);
		
		splice = getEnvelope(s, x1, x2);

		if (splice.first < splice.second) {
			right = x2;
		}
		else {
			left = x1;
		}

	}

	// ������������ �������� ����� ��������
	argmin = (right + left) / 2;
	splice = getEnvelope(s, argmin, argmin);
	min = splice.first;

	//if (s->p2->x - argmin < eps || argmin - s->p1->x < eps) {
	//	return pair<double, double>{ INFINITY, INFINITY };
	//}
	// �������� �������, ��� ����� ������ ������� ��� ��������� ���������
	currDiscrep = min;
	return pair<double, double>{ argmin, min };
}