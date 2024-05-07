#pragma once
#include "BaseMethod.h"
#include <fstream>
#include <iostream>
using namespace std;
class ParabolMethod : public BaseMethod {
public:
	//�����������
	ParabolMethod(shared_ptr<TProblem> _t, int _mode = 1) : BaseMethod(_t) {
		mode = _mode;
		fileNamePoints = "resultPointsParabol.txt";
		fileOperationCharacteristics = "resultIterationInformationParabol.txt";
		generateFibonacci(50);
	};
	virtual ~ParabolMethod() {};
protected:
	int mode; // ����� ������ (��������-������������� ��� �� �������)

	vector<double> Fibonacci;

	void setMode(int _mode) {
		mode = _mode;
	}

	void Test(const shared_ptr<Segment>& s);

	//��������� �������������� ���������
	virtual void calculateR(const shared_ptr<Segment>& s) override;

	// ��������-������������� �����
	void calculateRAnalytics(const shared_ptr<Segment>& s);

	// ����� ��������� �� �������
	void calculateRSections(const shared_ptr<Segment>& s);

	//��������� ��������� ������� �����������
	void calculateLowerLimit(const shared_ptr<Segment>& s) override;
	//��������� m1, m2, m3, m4 - ����� ������������
	double findLLoc(const shared_ptr<Segment>& s, const int i) const override;
private:
	//�������� ������� ����� ���������
	void generateFibonacci(int n = 100);

	// ������� ��������� ������ ������ �����������
	pair<double, double> getEnvelope(const shared_ptr<Segment>& s, double x1, double x2);

	//����� ���������
	pair<double, double> fibonacciMethod(const shared_ptr<Segment>& s, pair<double, double> seg);

	//����� ����� ����������� ���������
	vector<double> findRoots(const double a, const double b, const double c) const;

	//����� ������ ����������� ���������, ��� ����� ������ �������
	vector<pair<double, double>> getNegAreaOfConstraints(const shared_ptr<Segment>& s) const;

	//����� ������ ����������� ��������� � ������ ���������� ���������� �����
	vector<pair<double, double>> getNegAreaOfConstraints(const shared_ptr<Segment>& s, double param) const;

	//����� �������� �� ������� �������� ������������ (��� �������� - ����� ����)
	vector<pair<double, double>> validSetOfNegPar(const double a, const double b, const vector<double> roots) const;

	//����� �������� �� ������� �������� ������������ (��� �������� - ����� �����)
	vector<pair<double, double>> validSetOfPosPar(const double a, const double b, const vector<double> roots) const;

	//����� ����������� ����������� �������� v1, v2, v3
	vector<pair<double, double>> combine(const vector<pair<double, double>> v1, const vector<pair<double, double>> v2, const vector<pair<double, double>> v3) const;

	//����� ����������� �������� ������������� �������� �� ��������� a, b
	pair<double, double> minPos(const double x, const double y, const double a,
		                        const double b, const double c, const double arg_min) const;

	//����� ����������� �������� ������������� �������� �� ��������� a, b
	pair<double, double> minNeg(const double x, const double y, const double a, const double b, const double c) const;

	//����������� �������� ���������� 
	vector<pair<double, double>> intersectUnions(const vector<pair<double, double>> x, const vector<pair<double, double>> y) const;

	//����� ����������� �������� ����������� ������� ������� �� ��������� p
	pair<double, double> findMin(const shared_ptr<Segment> s, const pair<double, double> p) const;

	// TODO remove
	virtual void printViewData() override;
	//���������� ������� ��������� �������� �����������  
	double calculateMinorant(const shared_ptr<Segment>& s, double x);
};

