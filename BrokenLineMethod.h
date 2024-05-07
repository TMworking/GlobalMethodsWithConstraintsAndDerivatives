#pragma once
#include "BaseMethod.h"

class BrokenLineMethod : public BaseMethod {
public:
	//�����������
	BrokenLineMethod(shared_ptr<TProblem> _t) : BaseMethod(_t) {
		fileNamePoints = "resultPointsBrokenLine.txt";
		fileOperationCharacteristics = "resultIterationInformationBrokenLine.txt";
	};
	virtual ~BrokenLineMethod() {};
protected:
	//��������� �������������� ���������
	void calculateR(const shared_ptr<Segment>& s) override;

	//��������� ��������� ������� �����������
	void calculateLowerLimit(const shared_ptr<Segment>& s) override;

	double findLLoc(const shared_ptr<Segment>& s, const int i) const override;
};

