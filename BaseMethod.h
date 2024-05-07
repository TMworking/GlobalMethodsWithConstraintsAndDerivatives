#pragma once
#include <vector>
#include <iostream>
#include <limits>
#include <string>
#include <fstream>
#include <memory>
#include "Problems.h"
#define INFINITY std::numeric_limits<double>::infinity()
using namespace std;

//����� ���������
struct MyPoint {
	double x; //�������� ����������
	// ��� ������� ����������� ������������ � ������, � ������� �������� �������� ������� ��� �� �����������
	vector<double> f; //�������� �������
	vector<double> fd; //�������� �����������
	int numOfFuncs; //���������� ������� (������� + �����������)
	bool inArea; //����, ��������� �� ����� � ������� ���������� ��������

	//�����������
	MyPoint(double _x, int _numOfFuncs, shared_ptr<TProblem> t) : x(_x), numOfFuncs (_numOfFuncs) {
		f.reserve(numOfFuncs);
		fd.reserve(numOfFuncs);
		inArea = true;
		for (int i = 0; i < numOfFuncs + 1; i++) {
			f.push_back((*t)(x, i));
			fd.push_back((*t)(x, i, 'd'));
		}
		for (int i = 1; i < numOfFuncs + 1; i++) {
			if (f[i] > 0) {
				inArea = false;
				break;
			}
		}
	}
	//����������� �����������
	MyPoint(MyPoint& mp): x(mp.x), numOfFuncs(mp.numOfFuncs), inArea(mp.inArea) {
		f.reserve(numOfFuncs);
		fd.reserve(numOfFuncs);
		for (int i = 0; i < mp.numOfFuncs; ++i) {
			f[i] = mp.f[i];
			fd[i] = mp.fd[i];
		}
	}

	//����������� �����������, ������� �������� ����� ������������� �� mp.f � mp.fd � f � fd
	MyPoint(MyPoint&& mp) noexcept
		: numOfFuncs(mp.numOfFuncs), x(mp.x), f(mp.f), fd(mp.fd), inArea(mp.inArea) {
		mp.f.clear();
		mp.fd.clear();
	}

	//����������
	~MyPoint() {		
	}

	//�������� ������������
	MyPoint& operator= (const MyPoint& mp) {
		x = mp.x;
		numOfFuncs = mp.numOfFuncs;
		f.clear();
		fd.clear();
		f.reserve(numOfFuncs);
		fd.reserve(numOfFuncs);
		for (int i = 0; i < mp.numOfFuncs; ++i) {
			f[i] = mp.f[i];
			fd[i] = mp.fd[i];
		}
		return *this;
	}

	//�������� ������������ ������������, ������� �������� ����� ������������� �� mp.f � mp.fd � f � fd
	MyPoint& operator=(MyPoint&& mp) noexcept
	{
		// �������� �� ����������������
		if (&mp == this)
			return *this;

		// ������� ��, ��� � ����� ������� ����� ������� ��������� 
		f.clear();
		fd.clear();

		x = mp.x;
		numOfFuncs = mp.numOfFuncs;
		inArea = mp.inArea;
		// �������� ����� ������������� �� mp.f � mp.fd � f � fd
		f = move(mp.f);
		fd = move(mp.fd);
	
		return *this;
	}
};

//�����������
struct Segment {
	shared_ptr<MyPoint> p1; //������ ����� ���������
	shared_ptr<MyPoint> p2; //������ ����� ���������
	double length; //����� �������������
	double argmin; //����� x �������������� ���������
	double R; //���������-�������������� ���������
	//��������� ������ �����������
	vector<double> c;
	vector<double> w;
	vector<double> xt1; //x'
	vector<double> xt2; //x''
	vector<double> wa; 
	vector<double> ca;
	vector<double> wb;
	vector<double> cb;
	//������
	vector<double> L_loc; //��������� ������
	vector<double> L_care; //���������� ������
	vector<double> L; // ������ � ������� ����� � ����� ��������
	vector<double> L_glob; // ���������� ���������� ������
	int size;

	//�����������
	Segment(shared_ptr<MyPoint> _p1, shared_ptr<MyPoint> _p2, int _size) {
		R = INFINITY;
		argmin = 0.0;
		length = _p2->x - _p1->x;
		size = _size;
		p1 = _p1;
		p2 = _p2;
		L_loc.resize(size);
		L_care.resize(size);
		L.resize(size);
	    L_glob.resize(size);
		c.resize(size);
		w.resize(size);
		xt1.resize(size);
		xt2.resize(size);
		wa.resize(size);
		wb.resize(size);
		ca.resize(size);
		cb.resize(size);
	}

	//����������� �����������
	Segment(Segment& s) {
		p1 = s.p1;
		p2 = s.p2;
		size = s.size;
		R = s.R;
		argmin = s.argmin;
		length = s.length;
		L_loc.resize(size);
		L_care.resize(size);
		L.resize(size);
		L_glob.resize(size);
		c.resize(size);
		w.resize(size);
		xt1.resize(size);
		xt2.resize(size);
		wa.resize(size);
		wb.resize(size);
		ca.resize(size);
		cb.resize(size);
		for (int i = 0; i < size; i++) {
			L_loc[i] = s.L_loc[i];
			L_care[i] = s.L_care[i];
			L[i] = s.L[i];
			L_glob[i] = s.L_glob[i];
			c[i] = s.c[i];
			w[i] = s.w[i];
			xt2[i] = s.xt2[i];
			xt1[i] = s.xt1[i];
			wa[i] = s.wa[i];
			wb[i] = s.wb[i];
			ca[i] = s.ca[i];
			cb[i] = s.cb[i];
		}
	}

	//��������� ����������
	~Segment() = default;

	//������ � ����
	void print(string testInfoFile, int iteration) {
		std::ofstream testStream;
		testStream.open(testInfoFile, ios_base::app);

		//iteration
		testStream << "iteration=" << iteration << endl;
		testStream.precision(15);

		// left point
		testStream << " x1=" << p1->x;
		for (int i = 0; i < p1->numOfFuncs; i++) {
			testStream << " f[" << i << "]=" << p1->f[i];
			testStream << " f'[" << i << "]=" << p1->fd[i];
		}
		testStream << endl;

		//right point
		testStream << " x2=" << p2->x;
		for (int i = 0; i < p2->numOfFuncs; i++) {
			testStream << " f[" << i << "]=" << p2->f[i];
			testStream << " f'[" << i << "]=" << p2->fd[i];
		}
		testStream << endl;

		//characteristics
		testStream << " R=" << R;
		//testStream << " fMinGlob=" << fMinGlob;
		testStream << " argmin=" << argmin << endl;
		for (int i = 0; i < size; i++) {
			{
				testStream << " L_loc[" << i << "]=" << L_loc[i];
				testStream << " L_care[" << i << "]=" << L_care[i];
				testStream << " L_glob[" << i << "]=" << L_glob[i];

				testStream << " w[" << i << "]=" << w[i];
				testStream << " c[" << i << "]=" << c[i];
				testStream << " wa[" << i << "]=" << wa[i];
				testStream << " ca[" << i << "]=" << ca[i];
				testStream << " wb[" << i << "]=" << wb[i];
				testStream << " cb[" << i << "]=" << cb[i];
				testStream << " x'[" << i << "]=" << xt1[i];
				testStream << " x''[" << i << "]=" << xt2[i];
				testStream << endl;
			}
		}
		testStream << endl;
		testStream.close();
	}
};

// ������ ���������
enum class StopCriterion {
	CutLength,
	GeneratorClose
};

//��������� ���������� �����
struct Gamma {
	double g;
	double g_loc;
	double g_constraints;
	double g_loc_constraints;
	Gamma() :
		g(4.0f),
		g_loc(1.05f),
		g_constraints(3.0f),
		g_loc_constraints(1.05f) {}

	Gamma(double _g, double _g_loc,
		double _g_constraints, double _g_loc_constraints) :
		g(_g),
		g_loc(_g_loc),
		g_constraints(_g_constraints),
		g_loc_constraints(_g_loc_constraints) {}

	Gamma& operator=(const Gamma& gamma) {
		g = gamma.g;
		g_loc = gamma.g_loc;
		g_constraints = gamma.g_constraints;
		g_loc_constraints = gamma.g_loc_constraints;
		return *this;
	}
};

// ������� ����� �� ������� ���������� ����� �������� � � �������� �������� � ���� ������ �������
// ���������� ��������� �� ������ TProblem
class BaseMethod {
public:
	BaseMethod(shared_ptr<TProblem> _t);
	virtual ~BaseMethod();

	//���������� ���� ���������� �� ������ � ����
	void setLogNeeded(const bool _logNeeded);

	//���������� �������� ����������
	void setEps(const double _eps);

	//���������� ������������ ���������� ��������
	void setMaxIterationsCount(const int iterationsCount);

	//���������� �������� ���������
	void setStopCriterion(const StopCriterion stopCriterion);

	//���������� ���������
	void setParams(const int N, const Gamma gamma, const double v, const double d);

	void setNeedChangeGammaValue(const int step, const Gamma gamma, const int changeGam_Step)
	{
		changeGammaStep = step;
		gamma2 = gamma;
		needChangeGamma = true;
		changeGammaStep = changeGam_Step;
	}

	//�������� ���������� ��������
	int getIterationsCount() const;

	//��������
	void calculate();
	// TODO remove
	virtual void printViewData() {}; //��� ������ �������� ������� � �� ������ ������
	                            //��� ����������� ������������ �� ��������

protected:
	bool logNeeded = false; //������ � ����
	shared_ptr<TProblem> t = NULL;
	size_t constraintsCount; // ���������� �����������
	Gamma gamma1; //��������� ��� ���������� ������
	bool needChangeGamma = false; // ���������� �������������
	Gamma gamma2; // ��������� ��� ��������������� ���������
	int changeGammaStep = 5; // �������� ������������
	int N = 50; // ����� ����������� ��������� ����� * �������� �� + ��� ��������� ���������� ������
	double delta = 1.0e-10;
	double v = 0.2;
	double d = 0.25f; //���� �� ����� ���. ���������
	double eps = 1e-6;
	int iterationsCount = 0;
	int maxIterationCount = 500;
	StopCriterion stopCriterion = StopCriterion::GeneratorClose;

	vector<shared_ptr<Segment>> segments; // faster then std::list, the same as std::deque, � ������� �������� �� std::deque
	vector<shared_ptr<MyPoint>> points; // ������ �����, � ������� ����������� ���������, ���� �������� �� std::deque
	shared_ptr<Segment> st = NULL; //������� � ������� ���������� ���������� �� ������� ����
	vector<double> degree_of_increase;
	vector<double> L_glob;
	bool L_glob_changed = false;
	double fMinGlob; //f*- ������� ����� ������������ �������� �������
	bool newfMinGlob;//���� ��������� ��������� ���������� ��������

	double minDiscrep = INFINITY; // ���. ��������� �������� �������
	double currDiscrep = INFINITY; // �������� �������, ������� ���������� ������
	bool accVal = false; // ���� �� ���������� �������� G(X)<0

	string fileNamePoints; //����� � ���� ���� ����� ���������
	string fileOperationCharacteristics; //����� � ���� ��������� ���������� � ������ ��������

	// ���������� ���� �� ���������� ��������
	void defineAccVal(const shared_ptr<Segment>& s, const shared_ptr<MyPoint>& p);

	virtual bool isFullRecalcNeeded();
	//��������� ����� ����� � ����� ������� �� ���

	virtual void splitSegmentWithNewPoint(const shared_ptr<Segment>& s, shared_ptr<Segment>& s1, shared_ptr<Segment>& s2);
	virtual void recalcSegmentCharacteristics(const shared_ptr<Segment>& s);
	// ���������� ����� �� �������� ������������� ��� ���� ���������
	
	
	//��������� ��������� ������� ������� �����������, �����������, ������ ���� ����������� � ������� ��������
	virtual void calculateLowerLimit(const shared_ptr<Segment>& s) = 0;
	//��������� �������������� ���������, �����������, ������ ���� ����������� � ������� ��������
	virtual void calculateR(const shared_ptr<Segment>& s) = 0;
	virtual double findLLoc(const shared_ptr<Segment>& s, const int i) const = 0;

	//����������� ����������
	pair<double, double> intersectSegments(const pair<double, double> x, const pair<double, double> y) const;
	//���������� ������ ��������� �������
	void calculateL_new(const shared_ptr<Segment>& s);
	void calculateL_care(const shared_ptr<Segment>& s); // ��� ���������� ������� ������
	void calculateL_glob(const shared_ptr<Segment>& s, const int i, const double g, const double g_c);
	//��������� f*, �� ���� fMinGlob
	virtual void calculateFminGlob(shared_ptr<MyPoint> newPoint);
private:
	//�������� �������� ���������
	bool needStop() const;
	//�������� ����� �����
	void insertNewPoint(const shared_ptr<Segment>& s);
	void defineNextIteration();
	void recalcAdjacentSegments(const shared_ptr<Segment>& s);
	void calculateNewSegment(const shared_ptr<Segment>& s);
	double findLLocMax(const int i) const;
};
