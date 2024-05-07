#pragma once
#include "ParabolMethod.h"
#include <fstream>
#include <iostream>
using namespace std;

class TransformMethod : public ParabolMethod {

public:
	TransformMethod(shared_ptr<TProblem> _t) : ParabolMethod(_t) {
		fMinGlob = INFINITY;
		newfMinGlob = false;
		fileNamePoints = "resultPointsTransform.txt";
		fileOperationCharacteristics = "resultIterationInformationTransform.txt";
	};
	virtual ~TransformMethod() {};
protected:
	//��������� ��������������-R ���������
	void calculateR(const shared_ptr<Segment>& s) override;
	//������� ������� �� ��� �����
	void splitSegmentWithNewPoint(const shared_ptr<Segment>& s, shared_ptr<Segment>& s1, shared_ptr<Segment>& s2) override {
		BaseMethod::splitSegmentWithNewPoint(s, s1, s2);
		BaseMethod::calculateFminGlob(s1->p2);
	}
	//����� �� ������ ��������
	bool isFullRecalcNeeded() override {
		bool recalcNeeded = newfMinGlob || BaseMethod::isFullRecalcNeeded();
		newfMinGlob = false;
		return recalcNeeded;
	}
private:
	// TODO remove
	virtual void printViewData() override;
	//���������� ������� ��������� �������� ����������� 
	//  � ������������������ ������� �������: 
	double calculateMinorant(const shared_ptr<Segment>& s, double x);
};
